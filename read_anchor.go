package meta

import (
	"github.com/mingzhi/ncbiutils"
)

type MappedRead struct {
	Seq    []byte // mapped read sequence
	Pos    int    // start position
	Genome string // gneome name
}

// ReadAnchor maps reads to ortholog alignments.
// readMap: {GenomeName: SamRecords}
func ReadAnchor(readMap map[string]SamRecords, alignments []ncbiutils.SeqRecords) [][]MappedRead {
	groups := [][]MappedRead{}
	couldnot := make(map[string]bool)
	for _, aln := range alignments {
		mappedReads := []MappedRead{}
		for _, sr := range aln {
			genome := sr.Genome
			reads, found := readMap[genome]
			if found {
				start := reads.Search(sr.Loc.From)
				for {
					i := start
					if i >= len(reads) {
						break
					}

					r := reads[i]
					if r.Pos >= sr.Loc.To {
						break
					}

					s := map2Ref(r)
					mappedReads = append(mappedReads,
						MappedRead{Seq: s, Pos: r.Pos, Genome: genome})
					start++
				}
			} else {
				couldnot[sr.Genome] = true
			}
		}
		if len(mappedReads) > 3 {
			groups = append(groups, mappedReads)
		}

	}

	Info.Println(couldnot)

	return groups
}
