package meta

import (
	"math"
)

// Calculate total covariance for all mapped reads in a genome.
// records: SamRecords;
// genome: Genome;
// maxl: max length of correlations.
func CovGenome(records SamRecords, genome Genome, maxl int) *CovCalculator {
	cc := NewCovCalculator(maxl, true)
	for _, rec := range records {
		if findRefAcc(rec.Ref.Name()) == findRefAcc(genome.Accession) {
			read := map2Ref(rec) // mapped read to the reference genome.
			if rec.Pos+len(read) <= len(genome.Seq) {
				nucl := genome.Seq[rec.Pos : rec.Pos+len(read)]
				pos := genome.PosProfile[rec.Pos : rec.Pos+len(read)]
				subs := SubProfile(read, nucl, pos)

				for l := 0; l < maxl; l++ {
					for i := 0; i < len(subs)-l; i++ {
						x, y := subs[i], subs[i+l]
						if !math.IsNaN(x) && !math.IsNaN(y) {
							cc.Increment(l, x, y)
						}
					}
				}
			} else {
				Warn.Printf("%d, %d, %d\n", rec.Pos-1, rec.Pos+len(read), len(genome.PosProfile))
			}

		}
	}

	return cc
}

// Generate substitution profile according to the position profile.
func SubProfile(read, nucl, pos []byte) []float64 {
	subs := make([]float64, len(pos))
	for i, p := range pos {
		if p == FourFold {
			if read[i] == nucl[i] {
				subs[i] = 0
			} else {
				subs[i] = 1
			}
		} else {
			subs[i] = math.NaN()
		}
	}
	return subs
}
