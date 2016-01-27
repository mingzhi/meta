package main

import (
	"bytes"
	"flag"
	"fmt"
	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
	"github.com/mingzhi/biogo/feat/gff"
	"github.com/mingzhi/biogo/seq"
	"github.com/mingzhi/gomath/stat/correlation"
	"github.com/mingzhi/ncbiftp/genomes/profiling"
	"github.com/mingzhi/ncbiftp/taxonomy"
	"io"
	"log"
	"math"
	"os"
	"runtime"
)

// MappedRead contains the section of a read mapped to a reference genome.
type MappedRead struct {
	Pos  int
	Seq  []byte
	Qual []byte
}

type SubProfile struct {
	Pos     int
	Profile []float64
}

func (m MappedRead) Len() int {
	return len(m.Seq)
}

var MINBQ int
var MINMQ int

func main() {
	// Command variables.
	var bamFile string      // bam or sam file
	var genomeFile string   // genome file
	var gffFile string      // gff file
	var outFile string      // output file
	var maxl int            // max length of correlation
	var pos int             // position for calculation
	var codonTableID string // codon table ID
	var ncpu int            // number of CPUs
	// Parse command arguments.
	flag.IntVar(&maxl, "maxl", 100, "max length of correlations")
	flag.IntVar(&pos, "pos", 4, "position")
	flag.StringVar(&codonTableID, "codon", "11", "codon table ID")
	flag.IntVar(&ncpu, "ncpu", runtime.NumCPU(), "number of CPU for using")
	flag.IntVar(&MINBQ, "min-bq", 13, "min base quality")
	flag.IntVar(&MINMQ, "min-mq", 0, "min map quality")
	flag.Parse()
	// Print usage if the number of arguments is not satisfied.
	if flag.NArg() < 4 {
		log.Fatalln("Usage: go run calc_cr.go <pi file> <genome file> <gff file> <out file>")
	}
	bamFile = flag.Arg(0)
	genomeFile = flag.Arg(1)
	gffFile = flag.Arg(2)
	outFile = flag.Arg(3)
	runtime.GOMAXPROCS(ncpu)

	// Profile genome.
	// We need:
	// 1. genome file;
	// 2. gene features;
	// 3. condon table to identify four-fold degenerate sites.
	genome := readGenome(genomeFile)
	gffs := readGff(gffFile)
	codonTable := taxonomy.GeneticCodes()[codonTableID]
	profile := profiling.ProfileGenome(genome, gffs, codonTable)

	// Read sequence reads.
	_, readChan := readBamFile(bamFile)
	subProfileChan := slideReads(readChan)
	posType := convertPosType(pos)
	covs := calc(subProfileChan, profile, posType, maxl)
	write(covs, outFile)
}

// slideReads
func slideReads(readChan chan *sam.Record) chan SubProfile {
	subProfileChan := make(chan SubProfile)

	mappedReadArrChan := make(chan []MappedRead)
	go func() {
		defer close(mappedReadArrChan)
		mappedReadArr := []MappedRead{}
		for r := range readChan {
			if int(r.MapQ) > MINMQ {
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
			}
		}
	}()

	ncpu := runtime.GOMAXPROCS(0)
	done := make(chan bool)
	for i := 0; i < ncpu; i++ {
		go func() {
			for mappedReadArr := range mappedReadArrChan {
				a := mappedReadArr[0]
				mappedReadArr = mappedReadArr[1:]
				for _, b := range mappedReadArr {
					if b.Pos > a.Len()+a.Pos {
						break
					}
					subProfile := compareMappedReads(a, b)
					subProfileChan <- subProfile
				}
			}
			done <- true
		}()
	}

	go func() {
		defer close(subProfileChan)
		for i := 0; i < ncpu; i++ {
			<-done
		}
	}()

	return subProfileChan
}

// compareMappedReads compares two MappedReads in their overlapped part,
// and return a subsitution profile.
func compareMappedReads(a, b MappedRead) SubProfile {
	var subs []float64
	lag := b.Pos - a.Pos
	for j := 0; j < a.Len()-lag && j < b.Len(); j++ {
		i := j + lag
		d := math.NaN()
		if isATGC(a.Seq[i]) && isATGC(b.Seq[j]) {
			if int(a.Qual[i]) > MINBQ && int(b.Qual[j]) > MINBQ {
				if a.Seq[i] != b.Seq[j] {
					d = 1.0
				} else {
					d = 0.0
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

// calc
func calc(subProfileChan chan SubProfile, profile []profiling.Pos, posType byte, maxl int) (covs []*correlation.BivariateCovariance) {
	ncpu := runtime.GOMAXPROCS(0)
	covsChan := make(chan []*correlation.BivariateCovariance)
	for i := 0; i < ncpu; i++ {
		go func() {
			covs := []*correlation.BivariateCovariance{}
			for i := 0; i < maxl; i++ {
				covs = append(covs, correlation.NewBivariateCovariance(false))
			}

			for subProfile := range subProfileChan {
				for i := 0; i < len(subProfile.Profile); i++ {
					pos1 := subProfile.Pos + i
					x := subProfile.Profile[i]
					if checkPosType(posType, profile[pos1].Type) && !math.IsNaN(x) {
						for j := i; j < len(subProfile.Profile); j++ {
							pos2 := subProfile.Pos + j
							l := pos2 - pos1
							if l >= len(covs) {
								break
							} else {
								y := subProfile.Profile[j]
								if checkPosType(posType, profile[pos2].Type) && !math.IsNaN(y) {
									covs[l].Increment(x, y)
								}
							}

						}
					}

				}
			}
			covsChan <- covs
		}()
	}

	for i := 0; i < ncpu; i++ {
		covs1 := <-covsChan
		if len(covs) == 0 {
			covs = covs1
		} else {
			for i := 0; i < maxl; i++ {
				covs[i].Append(covs1[i])
			}
		}
	}

	return
}

// write
func write(covs []*correlation.BivariateCovariance, filename string) {
	w, err := os.Create(filename)
	if err != nil {
		log.Fatal(err)
	}
	defer w.Close()

	for i := 0; i < len(covs); i++ {
		v := covs[i].GetResult()
		n := covs[i].GetN()
		if n > 0 && !math.IsNaN(v) {
			w.WriteString(fmt.Sprintf("%d\t%g\t%d\n", i, v, n))
		}
	}
}

// readGenome read the genome file
// and return a genomic sequence.
func readGenome(filename string) []byte {
	f, err := os.Open(filename)
	if err != nil {
		panic(err)
	}
	defer f.Close()

	rd := seq.NewFastaReader(f)
	ss, err := rd.ReadAll()
	if err != nil {
		panic(err)
	}

	return ss[0].Seq
}

func readGff(filename string) []*gff.Record {
	f, err := os.Open(filename)
	if err != nil {
		panic(err)
	}
	defer f.Close()

	rd := gff.NewReader(f)
	ss, err := rd.ReadAll()
	if err != nil {
		panic(err)
	}

	records := []*gff.Record{}
	for _, s := range ss {
		if s.Feature == "CDS" {
			records = append(records, s)
		}
	}

	return records
}

type SamReader interface {
	Header() *sam.Header
	Read() (*sam.Record, error)
}

// ReadBamFile reads bam file, and return the header and a channel of sam records.
func readBamFile(fileName string) (h *sam.Header, c chan *sam.Record) {
	// Initialize the channel of sam records.
	c = make(chan *sam.Record)

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

		// Read and assign header.
		h = reader.Header()

		// Read sam records and send them to the channel,
		// until it hit an error, which raises a panic
		// if it is not a IO EOF.
		for {
			rec, err := reader.Read()
			if err != nil {
				if err != io.EOF {
					panic(err)
				}
				break
			}
			c <- rec
		}
		log.Println("Finished reading bam file!")
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

func checkPosType(posType, t1 byte) bool {
	isFirstPos := t1 == profiling.FirstPos
	isSecondPos := t1 == profiling.SecondPos
	isThirdPos := t1 == profiling.ThirdPos
	isFourFold := t1 == profiling.FourFold

	if posType == profiling.Coding {
		if isFirstPos || isSecondPos || isThirdPos || isFourFold {
			return true
		}
		return false
	}

	if posType == profiling.ThirdPos {
		if isThirdPos || isFourFold {
			return true
		}
		return false
	}

	return posType == t1
}

func convertPosType(pos int) byte {
	var p byte
	switch pos {
	case 0:
		p = profiling.NonCoding
	case 1:
		p = profiling.FirstPos
		break
	case 2:
		p = profiling.SecondPos
		break
	case 3:
		p = profiling.ThirdPos
		break
	case 4:
		p = profiling.FourFold
		break
	default:
		p = profiling.Coding
	}

	return p
}
