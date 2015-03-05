package meta

import (
	"code.google.com/p/biogo.bam/sam"
	"math"
	"runtime"
	"sort"
)

// Calculate total covariance for all mapped reads in a genome.
// records: SamRecords;
// genome: Genome;
// maxl: max length of correlations;
// pos: positions to be calculated.
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

// Calculate total covariance amoung all mapped reads.
// records: SamRecords;
// genome: Genome;
// maxl: max length of correlations;
// pos: postions to be calculated.
func CovReads(records SamRecords, genome Genome, maxl, pos int) (kc *KsCalculator, cc *CovCalculator) {
	// Find all records belongs to the genome.
	maxReadLength := 0
	founds := SamRecords{}
	for _, r := range records {
		if findRefAcc(r.Ref.Name()) == findRefAcc(genome.Accession) {
			founds = append(founds, r)
			if maxReadLength < r.Len() {
				maxReadLength = r.Len()
			}
		}
	}
	// Sort the records by left coordinate.
	sort.Sort(ByLeftCoordinate{founds})

	ncpu := runtime.GOMAXPROCS(0)
	// Create job channel.
	type job struct {
		i int
		r *sam.Record
	}
	jobs := make(chan job)
	go func() {
		for i, r := range founds {
			jobs <- job{i, r}
		}
		close(jobs)
	}()

	type result struct {
		cc *CovCalculator
		kc *KsCalculator
	}
	results := make(chan result)
	for i := 0; i < ncpu; i++ {
		go func() {
			// Prepare calculators.
			bias := true
			cc := NewCovCalculator(maxl, bias)
			kc := NewKsCalculator()
			for i, r1 := range founds {
				read1 := map2Ref(r1)
				for j := i - 1; j >= 0; j-- {
					r2 := founds[j]
					if r2.Pos+maxReadLength < r1.Pos {
						break
					}
					read2 := map2Ref(r2)
					// Check if overlap.
					if r2.Pos+len(read2) > r1.Pos {
						// Determine overlap regions (in genome coordinate).
						start := maxInt(r1.Pos, r2.Pos)
						end := minInt(r1.Pos+len(read1), r2.Pos+len(read2))
						// Prepare profile and read sequences.
						profile := genome.PosProfile[start:end]
						nucl1 := read1[start-r1.Pos : end-r1.Pos]
						nucl2 := read2[start-r2.Pos : end-r2.Pos]

						// nucl := genome.Seq[start:end]
						// Info.Printf("%d\t%d\t%d\n", len(nucl1), len(nucl2), len(profile))
						// Info.Printf("\n%s\n%s\n%s\n\n", string(nucl1), string(nucl2), string(nucl))

						// Subsitution profiling.
						subs := SubProfile(nucl1, nucl2, profile, pos)

						// Update cc and kc.
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
					}

				}
			}
			results <- result{cc, kc}
		}()
	}

	for i := 0; i < ncpu; i++ {
		res := <-results
		if i == 0 {
			cc = res.cc
			kc = res.kc
		} else {
			cc.Append(res.cc)
			kc.Append(res.kc)
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
