package main

import (
	"bytes"
	"fmt"
	"github.com/biogo/hts/sam"
	"sort"
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

func (s *SNP) String() string {
	bases := []byte{}
	for i := 0; i < len(s.Bases); i++ {
		bases = append(bases, s.Bases[i].Base)
	}
	return fmt.Sprintf("Pos: %d, Bases: %s", s.Pos, string(bases))
}

func Pileup(input chan *sam.Record) (output chan *SNP) {
	output = make(chan *SNP)

	go func() {
		defer close(output)
		// Initialize SNP map.
		m := make(map[int]*SNP)
		// Iterater over the records,
		// and for each record, obtain the bases,
		// and group them into the corresponding SNP.
		for r := range input {
			if int(r.MapQ) >= minMQ {
				s, q := Map2Ref(r)
				for i := 0; i < len(s); i++ {
					if s[i] != '*' {
						b := Base{}
						b.Base = s[i]
						b.Pos = r.Pos + i + 1
						b.Qual = q[i]
						b.ReadId = r.Name
						// if strings.Contains(r.Flags.String(), "1") {
						// 	b.ReadId += "_read1"
						// } else {
						// 	b.ReadId += "_read2"
						// }

						snp, found := m[b.Pos]
						if !found {
							snp = &SNP{}
							snp.Pos = b.Pos
						}
						snp.Bases = append(snp.Bases, &b)
						m[b.Pos] = snp
					}
				}

				completedSNPs := clearSNPMap(m, r.Pos)
				for _, snp := range completedSNPs {
					output <- snp
				}
			}
		}
	}()

	return
}

func clearSNPMap(m map[int]*SNP, minPos int) (snps []*SNP) {
	indices := []int{}
	for i := range m {
		indices = append(indices, i)
	}
	sort.Ints(indices)
	for _, i := range indices {
		if i < minPos-100 {
			s := m[i]
			snps = append(snps, FilterSNP(s))
			delete(m, i)
		}
	}

	return
}

// FilterSNP filters low quality bases,
// and overlapped bases in the pair-end reads.
func FilterSNP(s *SNP) *SNP {
	m := make(map[string]*Base)
	for i := 0; i < len(s.Bases); i++ {
		b := s.Bases[i]
		if int(b.Qual) >= minBQ {
			b1, found := m[b.ReadId]
			if found {
				if b.Base != b1.Base {
					delete(m, b.ReadId)
				}
			} else {
				m[b.ReadId] = b
			}
		}
	}

	ATGC := []byte{'A', 'T', 'G', 'C'}
	bases := []*Base{}
	for _, b := range m {
		b.Base = toUpper(b.Base)
		if bytes.Contains(ATGC, []byte{b.Base}) {
			bases = append(bases, b)
		} else {
			fmt.Printf("%d\t%s\t%v\n", b.Pos, b.ReadId, b.Base)
		}
	}
	s.Bases = bases

	return s
}

func toUpper(b byte) byte {
	var b1 byte
	switch b {
	case 'a':
		b1 = 'A'
	case 't':
		b1 = 'T'
	case 'g':
		b1 = 'G'
	case 'c':
		b1 = 'C'
	default:
		b1 = b
	}

	return b1
}
