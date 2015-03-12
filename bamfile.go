package meta

// BAM file operations.

import (
	"code.google.com/p/biogo.hts/bam"
	"code.google.com/p/biogo.hts/sam"
	"io"
	"log"
	"os"
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

// Read SAM file and return its header and records.
// NOT explicitly sorted.
func ReadSamFile(fileName string) (header *sam.Header, records []*sam.Record) {
	// Open sam file.
	f, err := os.Open(fileName)
	if err != nil {
		log.Panic(err)
	}
	defer f.Close()

	// Create sam reader,
	// and read the reference genomes.
	reader, err := sam.NewReader(f)
	if err != nil {
		log.Panic(err)
	}
	header = reader.Header()

	for {
		r, err := reader.Read()
		if err != nil {
			if err != io.EOF {
				log.Panic(err)
			} else {
				break
			}
		}
		records = append(records, r)
	}
	return
}

// Read BAM file and return its header and records.
// NOT explicitly sorted.
func ReadBamFile(fileName string) (header *sam.Header, records []*sam.Record) {
	// Open bam file.
	f, err := os.Open(fileName)
	if err != nil {
		log.Panic(err)
	}
	defer f.Close()

	// Create bam reader,
	// and read the reference genomes.
	rd := 0 // ignore this now.
	reader, err := bam.NewReader(f, rd)
	if err != nil {
		log.Panic(err)
	}
	header = reader.Header()

	for {
		r, err := reader.Read()
		if err != nil {
			if err != io.EOF {
				log.Panic(err)
			}
			break
		}
		records = append(records, r)
	}

	return
}

// Separate SAM records for different reference genomes.
// Return a map of genome reference name to records.
func SeparateSamRecords(refs []*sam.Reference, records SamRecords) map[string]SamRecords {
	m := make(map[string]SamRecords)
	for _, ref := range refs {
		m[ref.Name()] = FindSorted(ref.Name(), records)
	}
	return m
}

// Find and sort SamRecords of a reference.
func FindSorted(ref string, records SamRecords) SamRecords {
	founds := SamRecords{}
	for _, r := range records {
		if r.Ref.Name() == ref {
			founds = append(founds, r)
		}
	}
	Info.Println(ref)

	sort.Sort(ByLeftCoordinate{founds})

	return founds
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

func GetPairedEndReads(records SamRecords) PairedEndReads {
	// Sort records by name.
	sort.Sort(ByReadName{records})

	matedReads := []PairedEndRead{}
	var name string
	for _, r := range records {
		// Check if it has a mat
		if r.MateRef != nil && r.Ref == r.MateRef {
			if name != r.Name {
				matedReads = append(matedReads, PairedEndRead{})
				name = r.Name
			}

			lastIndex := len(matedReads) - 1
			if r.Pos < r.MatePos {
				matedReads[lastIndex].ReadLeft = r
			} else {
				matedReads[lastIndex].ReadRight = r
			}

			matedReads[lastIndex].Name = r.Name
		}
	}

	// double check
	matedReads2 := []PairedEndRead{}
	for _, r := range matedReads {
		if r.ReadLeft != nil && r.ReadRight != nil {
			matedReads2 = append(matedReads2, r)
		}
	}

	return matedReads2
}
