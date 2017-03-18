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

// MINBQ min bq
var MINBQ int

// MINMQ min mapping quality
var MINMQ int

// ShowProgress show progress.
var ShowProgress bool

func main() {

	// Command variables.
	var bamFile string // bam or sam file
	var outFile string // output file
	var maxl int       // max length of correlation
	var ncpu int       // number of CPUs
	// Parse command arguments.
	app := kingpin.New("meta_p2", "Calculate mutation correlation from bacterial metagenomic sequence data")
	app.Version("v0.1")
	bamFileArg := app.Arg("bamfile", "bam file").Required().String()
	outFileArg := app.Arg("outfile", "out file").Required().String()
	maxlFlag := app.Flag("maxl", "max len of correlations").Default("100").Int()
	ncpuFlag := app.Flag("ncpu", "number of CPUs").Default("0").Int()
	minBQFlag := app.Flag("minbq", "min base quality").Default("13").Int()
	minMQFlag := app.Flag("minmq", "min mapping quality").Default("30").Int()
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
	MINBQ = *minBQFlag
	MINMQ = *minMQFlag
	ShowProgress = *progressFlag

	runtime.GOMAXPROCS(ncpu)

	// Read sequence reads.
	refs, recordsChan := readBamFile(bamFile)
	codeTable := taxonomy.GeneticCodes()["11"]

	done := make(chan bool)
	p2Chan := make(chan []P2)
	for i := 0; i < ncpu; i++ {
		go func() {
			for records := range recordsChan {
				readsChan := slideReads(records)
				profileChan := compare(readsChan, codeTable)
				p2 := calc(profileChan, maxl)
				p2Chan <- p2
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

	numJob := len(refs)
	log.Printf("Number of references: %d\n", numJob)
	meanVars := collect(p2Chan, maxl, numJob)
	write(meanVars, outFile)
}

// slideReads
func slideReads(records []*sam.Record) chan []MappedRead {
	mappedReadArrChan := make(chan []MappedRead)
	go func() {
		defer close(mappedReadArrChan)

		totalDiscards := 0
		totalUsed := 0
		mappedReadArr := []MappedRead{}
		for _, r := range records {
			if r.MapQ != '*' && int(r.MapQ) < MINMQ {
				totalDiscards++
				continue
			}
			current := MappedRead{}
			current.Pos = r.Pos
			current.Seq, current.Qual = Map2Ref(r)
			mappedReadArr = append(mappedReadArr, current)
			if len(mappedReadArr) > 0 {
				a := mappedReadArr[0]
				if a.Pos+a.Len() < current.Pos {
					mappedReadArrChan <- mappedReadArr
					mappedReadArr = mappedReadArr[1:]
				}
			}
			totalUsed++
		}
	}()

	return mappedReadArrChan
}

func compare(readsChan chan []MappedRead, codeTable *taxonomy.GeneticCode) chan SubProfile {
	resChan := make(chan SubProfile)
	go func() {
		defer close(resChan)
		for reads := range readsChan {
			a := reads[0]
			for j := 1; j < len(reads); j++ {
				b := reads[j]
				if b.Pos > a.Len()+a.Pos {
					break
				}
				profile := compareMappedReads(a, b, codeTable)
				resChan <- profile
			}
		}
	}()
	return resChan
}

// compareMappedReads compares two MappedReads in their overlapped part,
// and return a subsitution profile.
func compareMappedReads(a, b MappedRead, codeTable *taxonomy.GeneticCode) SubProfile {
	var subs []float64
	lag := b.Pos - a.Pos
	for j := 0; j < a.Len()-lag && j < b.Len(); j++ {
		i := j + lag
		d := math.NaN()
		pos := j + b.Pos
		if (pos+1)%3 == 0 && j > 1 {
			if isATGC(a.Seq[i]) && isATGC(b.Seq[j]) {
				if int(a.Qual[i]) >= MINBQ && int(b.Qual[j]) >= MINBQ {
					codonA := string(a.Seq[i-2 : i+1])
					codonB := string(b.Seq[j-2 : j+1])
					aaA := codeTable.Table[codonA]
					aaB := codeTable.Table[codonB]
					if aaA == aaB {
						if a.Seq[i] != b.Seq[j] {
							d = 1.0
						} else {
							d = 0.0
						}
					}
				}
			}
		}
		subs = append(subs, d)
	}
	return SubProfile{Pos: b.Pos, Profile: subs}
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

type P2 struct {
	Total float64
	Count int
}

// calc
func calc(subProfileChan chan SubProfile, maxl int) (p2 []P2) {
	p2 = make([]P2, maxl)

	for subProfile := range subProfileChan {
		for i := 0; i < len(subProfile.Profile); i++ {
			x := subProfile.Profile[i]
			if !math.IsNaN(x) {
				for j := i; j < len(subProfile.Profile) && j-i < len(p2); j++ {
					l := j - i
					y := subProfile.Profile[j]
					if !math.IsNaN(y) {
						p2[l].Total += x * y
						p2[l].Count++
					}
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
		w.WriteString(fmt.Sprintf("%d,%g,%g,%d,%s,all\n", i, m, v, n, t))
	}
}

// SamReader is an interface for sam or bam reader.
type SamReader interface {
	Header() *sam.Header
	Read() (*sam.Record, error)
}

// ReadBamFile reads bam file, and return the header and a channel of sam records.
func readBamFile(fileName string) (refs []*sam.Reference, c chan []*sam.Record) {
	// Initialize the channel of sam records.
	c = make(chan []*sam.Record)

	// Create a new go routine to read the records.
	go func() {
		// Close the record channel when finished.
		defer close(c)

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
		refs = header.Refs()

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
					c <- records
					records = []*sam.Record{}
				}
				currentRefID = rec.RefID()
			}
			records = append(records, rec)
		}
		if len(records) > 0 {
			c <- records
		}
	}()

	return nil, c
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
