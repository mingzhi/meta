package main

import (
	"github.com/biogo/hts/sam"
	"runtime"
)

// BASE structure
type Base struct {
	Pos    int    // position it mapped to the reference genome.
	Base   byte   // the base pair.
	Qual   byte   // the quality score.
	ReadId string // the read id.
}

// Create sort interface
type Bases []*Base

func (b Bases) Len() int      { return len(b) }
func (b Bases) Swap(i, j int) { b[i], b[j] = b[j], b[i] }

// Sort by the read ID.
type ByReadId struct{ Bases }

func (b ByReadId) Less(i, j int) bool { return b.Bases[i].ReadId < b.Bases[j].ReadId }

// SNP structure
type SNP struct {
	Pos   int     // the position where it mapped to the reference genome.
	Bases []*Base // the mapped base pairs.
}

func Pileup(input chan *sam.Record) (output chan *SNP) {
	output = make(chan *SNP)

	ncpu := runtime.GOMAXPROCS(0)
	done := make(chan bool)
	for i := 0; i < ncpu; i++ {
		go func() {
			// Initialize SNP map.
			m := make(map[int]*SNP)
			// Iterater over the records,
			// and for each record, obtain the bases,
			// and group them into the corresponding SNP.
			for r := range input {
				s, q := Map2Ref(r)
				for i := 0; i < len(s); i++ {
					if s[i] != '*' {
						b := Base{}
						b.Base = s[i]
						b.Pos = r.Pos + i
						b.Qual = q[i]
						b.ReadId = r.Name
						snp, found := m[b.Pos]
						if !found {
							snp = &SNP{}
							snp.Pos = b.Pos
						}
						snp.Bases = append(snp.Bases, &b)
						m[b.Pos] = snp
					}
				}

				if len(m) > r.Len() {
					completedSNPs := ClearSNPMap(m, r.Pos)
					for _, snp := range completedSNPs {
						output <- snp
					}
				}
			}
			done <- true
		}()
	}

	// Create a new thread to monitor the job processing,
	// Close the output when all workers are done.
	go func() {
		defer close(output)
		for i := 0; i < ncpu; i++ {
			<-done
		}
	}()

	return
}

func ClearSNPMap(m map[int]*SNP, minPos int) (snps []*SNP) {
	for i, s := range m {
		if i < minPos {
			snps = append(snps, s)
			delete(m, i)
		}
	}

	return
}
