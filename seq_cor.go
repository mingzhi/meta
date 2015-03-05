package meta

import (
	"math"
)

// Calculate total covariance for all mapped reads in a genome.
// records: SamRecords;
// genome: Genome;
// maxl: max length of correlations.
func CovGenome(records SamRecords, genome Genome, maxl, pos int) (kc *KsCalculator, cc *CovCalculator) {
	cc = NewCovCalculator(maxl, true)
	kc = NewKsCalculator()
	for _, rec := range records {
		if findRefAcc(rec.Ref.Name()) == findRefAcc(genome.Accession) {
			read := map2Ref(rec) // mapped read to the reference genome.
			if rec.Pos+len(read) <= len(genome.Seq) {
				nucl := genome.Seq[rec.Pos : rec.Pos+len(read)]
				profile := genome.PosProfile[rec.Pos : rec.Pos+len(read)]
				subs := SubProfile(read, nucl, profile, pos)

				for l := 0; l < maxl; l++ {
					for i := 0; i < len(subs)-l; i++ {
						x, y := subs[i], subs[i+l]
						if !math.IsNaN(x) && !math.IsNaN(y) {
							cc.Increment(l, x, y)
						}

						if l == 0 {
							if !math.IsNaN(x) {
								kc.Increment(x)
							}
						}
					}
				}
			} else {
				Warn.Printf("%d, %d, %d\n", rec.Pos-1, rec.Pos+len(read), len(genome.PosProfile))
			}

		}
	}

	return
}

// Generate substitution profile according to the position profile.
func SubProfile(read, nucl, profile []byte, pos int) []float64 {
	// determine the codon position.
	var cp byte
	switch pos {
	case 1:
		cp = FirstPos
	case 2:
		cp = SecondPos
	case 3:
		cp = ThirdPos
	case 4:
		cp = FourFold
	}
	subs := make([]float64, len(profile))
	for i, p := range profile {
		var match bool
		if cp == ThirdPos {
			match = p == ThirdPos || p == FourFold
		} else {
			match = p == cp
		}
		if match {
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
