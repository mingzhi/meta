package reads

import (
	"github.com/biogo/hts/sam"
	"sort"
)

// A container for a group of SAM records.
// It implements sort.Interface for sorting and searching.
type SamRecords []*sam.Record

func (sr SamRecords) Len() int      { return len(sr) }
func (sr SamRecords) Swap(i, j int) { sr[i], sr[j] = sr[j], sr[i] }

func (sr SamRecords) Search(pos int) (index int) {
	index = sort.Search(len(sr), func(i int) bool { return sr[i].Pos >= pos })
	return
}

// A wrapper for sorting SAM records by left cordinate.
type ByLeftCoordinate struct{ SamRecords }

func (b ByLeftCoordinate) Less(i, j int) bool {
	return b.SamRecords[i].Pos < b.SamRecords[j].Pos
}

// A wrapper for sorting SAM records by read name.
type ByReadName struct{ SamRecords }

func (b ByReadName) Less(i, j int) bool {
	return b.SamRecords[i].Name < b.SamRecords[j].Name
}

// Container for paired-end reads.
type PairedEndRead struct {
	Name      string      // read name.
	ReadLeft  *sam.Record // left read.
	ReadRight *sam.Record // right read.
}

type PairedEndReads []PairedEndRead

func (p PairedEndReads) Len() int {
	return len(p)
}

func (p PairedEndReads) Swap(i, j int) {
	p[i], p[j] = p[j], p[i]
}

type ByNamePairedEndReads struct{ PairedEndReads }

func (p ByNamePairedEndReads) Less(i, j int) bool {
	return p.PairedEndReads[i].Name < p.PairedEndReads[j].Name
}

type ByLeftCoordinatePairedEndReads struct{ PairedEndReads }

func (p ByLeftCoordinatePairedEndReads) Less(i, j int) bool {
	return p.PairedEndReads[i].ReadLeft.Pos < p.PairedEndReads[j].ReadLeft.Pos
}

type ByRightCoordinatePairedEndReads struct{ PairedEndReads }

func (p ByRightCoordinatePairedEndReads) Less(i, j int) bool {
	r1 := p.PairedEndReads[i]
	r2 := p.PairedEndReads[j]
	return r1.ReadRight.Pos+r1.ReadRight.Len() < r2.ReadRight.Pos+r2.ReadRight.Len()
}
