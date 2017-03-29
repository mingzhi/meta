package main

import (
	"bytes"
	"fmt"
	"log"
	"math"
	"os"
	"runtime"

	"encoding/json"

	"bufio"

	"io"
	"strings"

	"github.com/biogo/hts/sam"
	"github.com/mingzhi/biogo/seq"
	"github.com/mingzhi/gomath/stat/desc/meanvar"
	"github.com/mingzhi/ncbiftp/taxonomy"
	"gopkg.in/alecthomas/kingpin.v2"
	"gopkg.in/cheggaaa/pb.v1"
)

// MappedRead contains the section of a read mapped to a reference genome.
type MappedRead struct {
	Pos  int
	Seq  []byte
	Qual []byte
}

// SubProfile Substitution/mutation profile.
type SubProfile struct {
	Pos     int
	Profile []float64
}

// Len return the lenght of a sequence.
func (m MappedRead) Len() int {
	return len(m.Seq)
}

// ShowProgress show progress.
var ShowProgress bool

// MinBaseQuality min base quality
var MinBaseQuality int

// MinMapQuality min map quality
var MinMapQuality int

// MinAlleleDepth min allele depth.
var MinAlleleDepth int

func main() {
	// Command variables.
	var bamFile string      // bam or sam file
	var outFile string      // output file
	var maxl int            // max length of correlation
	var ncpu int            // number of CPUs
	var minDepth int        // min depth
	var minCoverage float64 // min coveage
	var gffFile string      // gff file
	var corrResFile string  // corr result file.
	var geneFile string     // gene file.

	// Parse command arguments.
	app := kingpin.New("meta_p2", "Calculate mutation correlation from bacterial metagenomic sequence data")
	app.Version("v0.1")
	bamFileArg := app.Arg("bamfile", "bam file").Required().String()
	outFileArg := app.Arg("outfile", "out file").Required().String()
	maxlFlag := app.Flag("maxl", "max len of correlations").Default("100").Int()
	ncpuFlag := app.Flag("ncpu", "number of CPUs").Default("0").Int()
	minDepthFlag := app.Flag("min-depth", "min depth").Default("10").Int()
	minCoverageFlag := app.Flag("min-coverage", "min coverage").Default("0.8").Float64()
	progressFlag := app.Flag("progress", "show progress").Default("false").Bool()
	gffFileFlag := app.Flag("gff-file", "gff file").Default("").String()
	minBaseQFlag := app.Flag("min-base-qual", "min base quality").Default("30").Int()
	minMapQFlag := app.Flag("min-map-qual", "min mapping quality").Default("30").Int()
	corrResFileFlag := app.Flag("corr-res-file", "corr result file").Default("").String()
	geneFileFlag := app.Flag("gene-file", "gene file").Default("").String()
	minAlleleDepthFlag := app.Flag("min-allele-depth", "min allele depth").Default("0").Int()
	kingpin.MustParse(app.Parse(os.Args[1:]))

	bamFile = *bamFileArg
	outFile = *outFileArg
	maxl = *maxlFlag
	if *ncpuFlag == 0 {
		ncpu = runtime.NumCPU()
	} else {
		ncpu = *ncpuFlag
	}
	ShowProgress = *progressFlag
	minDepth = *minDepthFlag
	minCoverage = *minCoverageFlag
	gffFile = *gffFileFlag
	MinBaseQuality = *minBaseQFlag
	MinMapQuality = *minMapQFlag
	corrResFile = *corrResFileFlag
	geneFile = *geneFileFlag
	MinAlleleDepth = *minAlleleDepthFlag

	runtime.GOMAXPROCS(ncpu)

	// Read sequence reads.
	var header *sam.Header
	var recordsChan chan GeneSamRecords
	if gffFile != "" {
		gffRecMap := readGffs(gffFile)
		header, recordsChan = readStrainBamFile(bamFile, gffRecMap)
	} else {
		header, recordsChan = readPanGenomeBamFile(bamFile)
	}

	var geneSet map[string]bool
	if geneFile != "" {
		genes := readLines(geneFile)
		for _, gene := range genes {
			geneSet[gene] = true
		}
	}

	codeTable := taxonomy.GeneticCodes()["11"]

	done := make(chan bool)
	p2Chan := make(chan CorrResults)
	for i := 0; i < ncpu; i++ {
		go func() {
			for geneRecords := range recordsChan {
				if geneFile != "" {
					if !geneSet[geneRecords.ID] {
						continue
					}
				}
				geneLen := geneRecords.End - geneRecords.Start
				gene := pileupCodons(geneRecords)
				ok := checkCoverage(gene, geneLen, minDepth, minCoverage)
				if ok {
					p2 := calcP2(gene, 21, minDepth, codeTable)
					p4 := calcP4(gene, maxl, minDepth, codeTable)
					p2 = append(p2, p4...)
					p2Chan <- CorrResults{Results: p2, GeneID: geneRecords.ID, GeneLenth: geneLen, ReadNum: len(geneRecords.Records)}
				}
			}
			done <- true
		}()
	}

	go func() {
		defer close(p2Chan)
		for i := 0; i < ncpu; i++ {
			<-done
		}
	}()

	var corrResEncoder *json.Encoder
	if corrResFile != "" {
		f, err := os.Create(corrResFile)
		if err != nil {
			log.Panic(err)
		}
		defer f.Close()
		corrResEncoder = json.NewEncoder(f)
	}
	collector := NewCollector()
	for corrResults := range p2Chan {
		collector.Add(corrResults)
		if corrResFile != "" {
			if err := corrResEncoder.Encode(corrResults); err != nil {
				log.Panic(err)
			}
		}
	}

	numJob := len(header.Refs())
	log.Printf("Number of references: %d\n", numJob)
	w, err := os.Create(outFile)
	if err != nil {
		panic(err)
	}
	defer w.Close()

	w.WriteString("l,m,v,n,t\n")
	results := collector.Results()
	for _, res := range results {
		w.WriteString(fmt.Sprintf("%d,%g,%g,%d,%s\n",
			res.Lag, res.Value, res.Variance, res.Count, res.Type))
	}
}

// pileupCodons pileup codons of a list of reads at a gene.
func pileupCodons(geneRecords GeneSamRecords) (codonGene *CodonGene) {
	codonGene = NewCodonGene()
	for _, read := range geneRecords.Records {
		if int(read.MapQ) < MinMapQuality || len(read.Cigar) > 1 || read.Cigar[0].Type() != sam.CigarMatch {
			continue
		}
		codonArray := getCodons(read, geneRecords.Start, geneRecords.Strand)
		for _, codon := range codonArray {
			codonGene.AddCodon(codon)
		}
	}

	return
}

// getCodons split a read into a list of Codon.
func getCodons(read *sam.Record, offset, strand int) (codonArray []Codon) {
	// get the mapped sequence of the read onto the reference.
	mappedSeq, _ := Map2Ref(read)
	for i := 2; i < len(mappedSeq); {
		if (read.Pos+i-offset+1)%3 == 0 {
			codonSeq := mappedSeq[i-2 : i+1]
			genePos := (read.Pos+i-offset+1)/3 - 1
			if genePos >= 0 {
				if strand == -1 {
					codonSeq = seq.Reverse(seq.Complement(codonSeq))
				}
				codon := Codon{ReadID: read.Name, Seq: string(codonSeq), GenePos: genePos}
				codonArray = append(codonArray, codon)
			}
			i += 3
		} else {
			i++
		}
	}

	return
}

func isATGC(b byte) bool {
	if b == 'A' {
		return true
	} else if b == 'T' {
		return true
	} else if b == 'C' {
		return true
	} else if b == 'G' {
		return true
	}

	return false
}

// P2 stores p2 calculation results.
type P2 struct {
	Total float64
	Count int
}

// doubleCount count codon pairs.
func doubleCount(nc *NuclCov, codonPairArray []CodonPair) {
	for _, cp := range codonPairArray {
		a := cp.A.Seq[2]
		b := cp.B.Seq[2]
		nc.Add(a, b)
	}
}

func calcP2(gene *CodonGene, maxl, minDepth int, codeTable *taxonomy.GeneticCode) (p2Res []CorrResult) {
	gene.SortCodonByReadID()
	alphabet := []byte{'A', 'T', 'G', 'C'}
	for i := 0; i < gene.Len(); i++ {
		for j := i; j < gene.Len(); j++ {
			codonPairRaw := gene.PairCodonAt(i, j)
			if len(codonPairRaw) < 2 {
				continue
			}
			lag := codonPairRaw[0].B.GenePos - codonPairRaw[0].A.GenePos
			if lag < 0 {
				lag = -lag
			}
			if lag >= maxl {
				break
			}

			splittedCodonPairs := SynoumousSplitCodonPairs(codonPairRaw, codeTable)
			for _, synPairs := range splittedCodonPairs {
				if len(synPairs) > minDepth {
					nc := NewNuclCov(alphabet)
					doubleCount(nc, synPairs)

					for len(p2Res) <= lag {
						p2Res = append(p2Res, CorrResult{Type: "P2", Lag: len(p2Res)})
					}
					xy, _, _, n := nc.Cov11(MinAlleleDepth)
					p2Res[lag].Count += int64(n)
					p2Res[lag].Value += xy
				}
			}
		}
	}

	return
}

func calcP4(gene *CodonGene, maxl, minDepth int, codeTable *taxonomy.GeneticCode) (p4Res []CorrResult) {
	gene.SortCodonByReadID()
	var valueArray []float64
	var countArray []int
	var posArray []int
	for i := 0; i < gene.Len(); i++ {
		value, count := autoCov(gene, i, minDepth, codeTable)
		if count > 0 {
			pos := gene.CodonPiles[i].Codons[0].GenePos
			valueArray = append(valueArray, value)
			countArray = append(countArray, count)
			posArray = append(posArray, pos)
		}

	}
	for i := 0; i < len(valueArray); i++ {
		value1 := valueArray[i]
		count1 := countArray[i]
		xbar := value1 / float64(count1)
		for j := i; j < len(valueArray); j++ {
			value2 := valueArray[j]
			count2 := countArray[j]
			ybar := value2 / float64(count2)
			lag := posArray[j] - posArray[i]
			if lag < 0 {
				lag = -lag
			}
			if lag >= maxl {
				break
			}
			for len(p4Res) <= lag {
				p4Res = append(p4Res, CorrResult{Type: "P4", Lag: len(p4Res)})
			}
			p4Res[lag].Value += xbar * ybar
			p4Res[lag].Count++
		}
	}

	return
}

func autoCov(gene *CodonGene, i, minDepth int, codeTable *taxonomy.GeneticCode) (value float64, count int) {
	alphabet := []byte{'A', 'T', 'G', 'C'}
	codonPairRaw := gene.PairCodonAt(i, i)
	if len(codonPairRaw) < 2 {
		return
	}
	lag := codonPairRaw[0].B.GenePos - codonPairRaw[0].A.GenePos
	if lag < 0 {
		lag = -lag
	}

	splittedCodonPairs := SynoumousSplitCodonPairs(codonPairRaw, codeTable)
	for _, synPairs := range splittedCodonPairs {
		if len(synPairs) > minDepth {
			nc := NewNuclCov(alphabet)
			doubleCount(nc, synPairs)

			xy, _, _, n := nc.Cov11(MinAlleleDepth)
			value += xy
			count += n
		}
	}
	return
}

// collect
func collect(p2Chan chan []P2, maxl, numJob int) (meanVars []*meanvar.MeanVar) {
	meanVars = []*meanvar.MeanVar{}
	for i := 0; i < maxl; i++ {
		meanVars = append(meanVars, meanvar.New())
	}

	var pbar *pb.ProgressBar
	if ShowProgress {
		pbar = pb.StartNew(numJob)
		defer pbar.Finish()
	}

	for p2 := range p2Chan {
		for i := range p2 {
			n := p2[i].Count
			v := p2[i].Total
			if n > 100 {
				meanVars[i].Increment(v / float64(n))
			}
		}
		if ShowProgress {
			pbar.Increment()
		}
	}

	return
}

// write
func write(meanVars []*meanvar.MeanVar, filename string) {
	w, err := os.Create(filename)
	if err != nil {
		log.Fatal(err)
	}
	defer w.Close()

	w.WriteString("l,m,v,n,t,b\n")
	ks := 0.0
	for i := 0; i < len(meanVars); i++ {
		m := meanVars[i].Mean.GetResult()
		v := meanVars[i].Var.GetResult()
		n := meanVars[i].Mean.GetN()
		if n == 0 || math.IsNaN(v) {
			continue
		}
		t := "P2"
		if i == 0 {
			t = "Ks"
			ks = m
		} else {
			m = m / ks
		}
		w.WriteString(fmt.Sprintf("%d,%g,%g,%d,%s,all\n", i*3, m, v, n, t))
	}
}

// Map2Ref Obtains a read mapping to the reference genome.
func Map2Ref(r *sam.Record) (s []byte, q []byte) {
	p := 0                 // position in the read sequence.
	read := r.Seq.Expand() // read sequence.
	qual := r.Qual
	for _, c := range r.Cigar {
		switch c.Type() {
		case sam.CigarMatch, sam.CigarMismatch, sam.CigarEqual:
			s = append(s, read[p:p+c.Len()]...)
			q = append(q, qual[p:p+c.Len()]...)
			p += c.Len()
		case sam.CigarInsertion, sam.CigarSoftClipped, sam.CigarHardClipped:
			p += c.Len()
		case sam.CigarDeletion, sam.CigarSkipped:
			for i := 0; i < c.Len(); i++ {
				s = append(s, '-')
				q = append(q, 0)
			}
		}
	}

	s = bytes.ToUpper(s)

	for i, a := range q {
		if int(a) < MinBaseQuality {
			s[i] = '-'
		}
	}

	return
}

func checkCoverage(gene *CodonGene, geneLen, minDepth int, minCoverage float64) (ok bool) {
	num := 0
	for _, pile := range gene.CodonPiles {
		if pile.Len() > minDepth {
			num++
		}
	}
	coverage := float64(num) / float64(geneLen/3)
	ok = coverage > minCoverage
	return
}

// readLines return all trimmed lines.
func readLines(filename string) []string {
	f, err := os.Open(filename)
	if err != nil {
		log.Panic(err)
	}
	defer f.Close()

	rd := bufio.NewReader(f)
	var lines []string
	for {
		line, err := rd.ReadString('\n')
		if err != nil {
			if err != io.EOF {
				log.Panic(err)
			}
			break
		}
		lines = append(lines, strings.TrimSpace(line))
	}
	return lines
}
