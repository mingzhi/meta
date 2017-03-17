package main

import (
	"bytes"
	"flag"
	"fmt"
	"io"
	"log"
	"math"
	"os"
	"runtime"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
	"github.com/mingzhi/gomath/stat/correlation"
	"github.com/mingzhi/gomath/stat/desc/meanvar"
	"github.com/mingzhi/ncbiftp/taxonomy"
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

func main() {
	// Command variables.
	var bamFile string      // bam or sam file
	var outFile string      // output file
	var maxl int            // max length of correlation
	var codonTableID string // codon table ID
	var ncpu int            // number of CPUs
	// Parse command arguments.
	flag.IntVar(&maxl, "maxl", 100, "max length of correlations")
	flag.StringVar(&codonTableID, "codon", "11", "codon table ID")
	flag.IntVar(&ncpu, "ncpu", runtime.NumCPU(), "number of CPU for using")
	flag.IntVar(&MINBQ, "min-bq", 13, "min base quality")
	flag.IntVar(&MINMQ, "min-mq", 0, "min map quality")
	flag.Parse()
	// Print usage if the number of arguments is not satisfied.
	if flag.NArg() < 2 {
		log.Fatalln("Usage: go run calc_cr.go <pi file> <genome file> <gff file> <out file>")
	}
	bamFile = flag.Arg(0)
	outFile = flag.Arg(1)
	runtime.GOMAXPROCS(ncpu)

	// Read sequence reads.
	recordsChan := readBamFile(bamFile)
	codeTable := taxonomy.GeneticCodes()["11"]

	done := make(chan bool)
	covsChan := make(chan []*correlation.BivariateCovariance)
	for i := 0; i < ncpu; i++ {
		go func() {
			for records := range recordsChan {
				readsChan := slideReads(records)
				profileChan := compare(readsChan, codeTable)
				covs := calc(profileChan, maxl)
				covsChan <- covs
			}
			done <- true
		}()
	}

	go func() {
		defer close(covsChan)
		for i := 0; i < ncpu; i++ {
			<-done
		}
	}()

	meanVars := collect(covsChan, maxl)
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
			if int(r.MapQ) > MINMQ && int(r.MapQ) < 51 {
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
			} else {
				totalDiscards++
			}
		}
		log.Printf("Total discard reads: %d\n", totalDiscards)
		log.Printf("Total used reads: %d\n", totalUsed)
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
				if int(a.Qual[i]) > MINBQ && int(b.Qual[j]) > MINBQ {
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
			subs = append(subs, d)
		}

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

// calc
func calc(subProfileChan chan SubProfile, maxl int) (covs []*correlation.BivariateCovariance) {

	covs = []*correlation.BivariateCovariance{}
	for i := 0; i < maxl; i++ {
		covs = append(covs, correlation.NewBivariateCovariance(false))
	}

	for subProfile := range subProfileChan {
		for i := 0; i < len(subProfile.Profile); i++ {
			pos1 := subProfile.Pos + i
			x := subProfile.Profile[i]
			if !math.IsNaN(x) {
				for j := i; j < len(subProfile.Profile); j++ {
					pos2 := subProfile.Pos + j
					l := pos2 - pos1
					if l >= len(covs) {
						break
					} else {
						y := subProfile.Profile[j]
						if !math.IsNaN(y) {
							covs[l].Increment(x, y)
						}
					}

				}
			}

		}
	}
	return

}

// collect
func collect(covsChan chan []*correlation.BivariateCovariance, maxl int) (meanVars []*meanvar.MeanVar) {
	meanVars = []*meanvar.MeanVar{}
	for i := 0; i < maxl; i++ {
		meanVars = append(meanVars, meanvar.New())
	}

	for covs := range covsChan {
		for i := range covs {
			c := covs[i]
			v := c.GetResult()
			if !math.IsNaN(v) {
				meanVars[i].Increment(v)
			}
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

	w.WriteString("l,m,v,n\n")
	for i := 0; i < len(meanVars); i++ {
		m := meanVars[i].Mean.GetResult()
		v := meanVars[i].Var.GetResult()
		n := meanVars[i].Mean.GetN()
		w.WriteString(fmt.Sprintf("%d,%g,%g,%d\n", i, m, v, n))
	}
}

// SamReader is an interface for sam or bam reader.
type SamReader interface {
	Header() *sam.Header
	Read() (*sam.Record, error)
}

// ReadBamFile reads bam file, and return the header and a channel of sam records.
func readBamFile(fileName string) chan []*sam.Record {
	// Initialize the channel of sam records.
	c := make(chan []*sam.Record)

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
		log.Println("Finished reading bam file!")
	}()

	return c
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
