package main

import (
	"bytes"
	"fmt"
	"io"
	"log"
	"math"
	"os"
	"runtime"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
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

func main() {
	// Command variables.
	var bamFile string // bam or sam file
	var outFile string // output file
	var maxl int       // max length of correlation
	var ncpu int       // number of CPUs
	var minDepth int   // min depth
	var minCoverage float64
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

	runtime.GOMAXPROCS(ncpu)

	// Read sequence reads.
	headerChan, recordsChan := readBamFile(bamFile)
	header := <-headerChan

	codeTable := taxonomy.GeneticCodes()["11"]

	done := make(chan bool)
	p2Chan := make(chan []CorrResult)
	for i := 0; i < ncpu; i++ {
		go func() {
			for records := range recordsChan {
				geneLen := records[0].Ref.Len()
				gene := pileupCodons(records, 0)
				ok := checkCoverage(gene, geneLen, minDepth, minCoverage)
				if ok {
					p2, p4 := calcP2(gene, maxl, minDepth, codeTable)
					p2Chan <- p2
					p2Chan <- p4
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

	collector := NewCollector()
	for corrResults := range p2Chan {
		collector.Add(CorrResults{Results: corrResults})
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
func pileupCodons(records []*sam.Record, offset int) (codonGene *CodonGene) {
	codonGene = NewCodonGene()
	for _, read := range records {
		codonArray := getCodons(read, offset)
		for _, codon := range codonArray {
			codonGene.AddCodon(codon)
		}
	}

	return
}

// getCodons split a read into a list of Codon.
func getCodons(read *sam.Record, offset int) (codonArray []Codon) {
	// get the mapped sequence of the read onto the reference.
	mappedSeq, _ := Map2Ref(read)
	for i := 2; i < len(mappedSeq); {
		if (read.Pos+i-offset+1)%3 == 0 {
			codonSeq := string(mappedSeq[i-2 : i+1])
			genePos := (read.Pos+i-offset+1)/3 - 1
			codon := Codon{ReadID: read.Name, Seq: codonSeq, GenePos: genePos}
			codonArray = append(codonArray, codon)
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

func calcP2(gene *CodonGene, maxl, minDepth int, codeTable *taxonomy.GeneticCode) (p2Res, p4Res []CorrResult) {
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
						p4Res = append(p4Res, CorrResult{Type: "P4", Lag: len(p4Res)})
					}
					xy, x, y, n := nc.Cov11()
					p2Res[lag].Count += n
					p2Res[lag].Value += xy

					xbar := float64(x) / float64(n)
					ybar := float64(y) / float64(n)
					p4Res[lag].Value += xbar * ybar
					p4Res[lag].Count++
				}
			}
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

// SamReader is an interface for sam or bam reader.
type SamReader interface {
	Header() *sam.Header
	Read() (*sam.Record, error)
}

// ReadBamFile reads bam file, and return the header and a channel of sam records.
func readBamFile(fileName string) (headerChan chan *sam.Header, recordChan chan []*sam.Record) {
	// Initialize the channel of sam records.
	headerChan = make(chan *sam.Header)
	recordChan = make(chan []*sam.Record)

	// Create a new go routine to read the records.
	go func() {
		// Close the record channel when finished.
		defer close(headerChan)
		defer close(recordChan)

		// Open file stream, and close it when finished.
		f, err := os.Open(fileName)
		if err != nil {
			panic(err)
		}
		defer f.Close()

		var reader SamReader
		if fileName[len(fileName)-3:] == "bam" {
			bamReader, err := bam.NewReader(f, 0)
			if err != nil {
				panic(err)
			}
			defer bamReader.Close()
			reader = bamReader
		} else {
			reader, err = sam.NewReader(f)
			if err != nil {
				panic(err)
			}
		}

		header := reader.Header()
		headerChan <- header

		// Read sam records and send them to the channel,
		// until it hit an error, which raises a panic
		// if it is not a IO EOF.
		currentRefID := -1
		var records []*sam.Record
		for {
			rec, err := reader.Read()
			if err != nil {
				if err != io.EOF {
					panic(err)
				}
				break
			}
			if currentRefID == -1 {
				currentRefID = rec.RefID()
			}
			if rec.RefID() != currentRefID {
				if len(records) > 0 {
					recordChan <- records
					records = []*sam.Record{}
				}
				currentRefID = rec.RefID()
			}
			records = append(records, rec)
		}
		if len(records) > 0 {
			recordChan <- records
		}
	}()

	return
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
				s = append(s, '*')
				q = append(q, 0)
			}
		}
	}

	s = bytes.ToUpper(s)

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
