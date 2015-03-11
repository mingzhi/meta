package meta

import (
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
	// Obtain paired-end reads.
	matedReads := GetPairedEndReads(records)
	// Info.Printf("%s, paired end reads: %d\n", genome.Accession, len(matedReads))
	sort.Sort(ByLeftCoordinatePairedEndReads{matedReads})
	// Prepare jobs.
	type job struct {
		r PairedEndRead
	}
	jobs := make(chan job)
	go func() {
		for _, r := range matedReads {
			sameRef := r.ReadLeft.Ref.Name() == r.ReadRight.Ref.Name()
			matchRef := FindRefAcc(r.ReadLeft.Ref.Name()) == FindRefAcc(genome.Accession)
			if sameRef && matchRef {
				jobs <- job{r}
			}
		}
		close(jobs)
	}()

	// Running jobs and send results to a chan.
	type result struct {
		cc *CovCalculator
		kc *KsCalculator
	}
	results := make(chan result)
	ncpu := runtime.GOMAXPROCS(0)
	for i := 0; i < ncpu; i++ {
		go func() {
			cc := NewCovCalculator(maxl, true)
			kc := NewKsCalculator()
			for j := range jobs {
				rec := j.r
				// mapped paired end read to the reference genome.
				read := mated2Ref(rec)
				if rec.ReadLeft.Pos+len(read) <= len(genome.Seq) {
					start := rec.ReadLeft.Pos
					end := rec.ReadLeft.Pos + len(read)
					nucl := genome.Seq[start:end]
					profile := genome.PosProfile[start:end]
					subs := SubProfile(read, nucl, profile, pos)
					// Info.Println("...")
					// Info.Printf("%v\n", rec.Name)
					// Info.Println(len(read))
					// Info.Printf("%d\t%d\n", rec.ReadLeft.Len(), rec.ReadRight.Len())
					// Info.Println(string(nucl))
					// Info.Println(string(read))
					// Info.Println("###")

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
					Warn.Printf("%d, %d, %d\n", rec.ReadLeft.Pos-1, rec.ReadLeft.Pos+len(read), len(genome.PosProfile))
				}
			}
			results <- result{cc, kc}
		}()
	}

	// Receive results from the chan.
	for i := 0; i < ncpu; i++ {
		r := <-results
		if i == 0 {
			cc = r.cc
			kc = r.kc
		} else {
			cc.Append(r.cc)
			kc.Append(r.kc)
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
	// Obtain paired-end reads.
	matedReads := GetPairedEndReads(records)
	founds := PairedEndReads{}
	maxReadLength := 0
	for _, r := range matedReads {
		if FindRefAcc(r.ReadLeft.Ref.Name()) == FindRefAcc(genome.Accession) {
			founds = append(founds, r)
			readLen := r.ReadRight.Pos + r.ReadRight.Len() - r.ReadLeft.Pos
			if maxReadLength < readLen {
				maxReadLength = readLen
			}
		}
	}
	// Info.Printf("%s, paired end reads: %d\n", genome.Accession, len(matedReads))
	sort.Sort(ByLeftCoordinatePairedEndReads{founds})

	ncpu := runtime.GOMAXPROCS(0)
	// Create job channel.
	type job struct {
		i int
		r PairedEndRead
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
				read1 := mated2Ref(r1)
				for j := i - 1; j >= 0; j-- {
					r2 := founds[j]
					if r2.ReadLeft.Pos+maxReadLength < r1.ReadLeft.Pos {
						break
					}
					read2 := mated2Ref(r2)
					// Check if overlap.
					if r2.ReadLeft.Pos+len(read2) > r1.ReadLeft.Pos {
						// Determine overlap regions (in genome coordinate).
						start := maxInt(r1.ReadLeft.Pos, r2.ReadLeft.Pos)
						end := minInt(r1.ReadLeft.Pos+len(read1), r2.ReadLeft.Pos+len(read2))
						// Prepare profile and read sequences.
						profile := genome.PosProfile[start:end]
						nucl1 := read1[start-r1.ReadLeft.Pos : end-r1.ReadLeft.Pos]
						nucl2 := read2[start-r2.ReadLeft.Pos : end-r2.ReadLeft.Pos]

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
		valid := read[i] != '*' && nucl[i] != '*'
		if match && valid {
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
