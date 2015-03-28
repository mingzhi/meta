package reads

import (
	"github.com/biogo/hts/sam"
	"sort"
)

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

	sort.Sort(ByLeftCoordinate{founds})

	return founds
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
