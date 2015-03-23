package meta

import (
	"github.com/mingzhi/ncbiutils"
	"math"
	"runtime"
)

type CovGenomesFunc func(alignments []ncbiutils.SeqRecords, genome Genome, maxl, pos int) (kc *KsCalculator, cc *CovCalculator)

// Calculate correlation of substituions in reads,
// by comparing them to the reference genome.
// records: SamRecords;
// genome: Genome;
// maxl: max length of correlations;
// pos: positions to be calculated.
func CovReadsGenome(matedReads PairedEndReads, genome Genome, maxl, pos int) (kc *KsCalculator, cc *CovCalculator) {
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
					SubCorr(subs, cc, kc, maxl)
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

// Calculate correlation of substituions in reads,
// by comparing reads to reads.
// records: SamRecords;
// genome: Genome;
// maxl: max length of correlations;
// pos: postions to be calculated.
func CovReadsReads(matedReads PairedEndReads, genome Genome, maxl, pos int) (kc *KsCalculator, cc *CovCalculator) {
	ncpu := runtime.GOMAXPROCS(0)
	// Create job channel.
	type job struct {
		r1, r2 PairedEndRead
	}
	jobs := make(chan job)
	go func() {

		for i := 0; i < len(matedReads); i++ {
			r1 := matedReads[i]
			for j := i - 1; j >= 0; j-- {
				r2 := matedReads[j]
				if r2.ReadRight.Pos+r2.ReadRight.Len() < r1.ReadLeft.Pos {
					break
				} else {
					jobs <- job{r1, r2}
				}
			}
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
			// prepare calculators.
			biasCorrection := true
			cc := NewCovCalculator(maxl, biasCorrection)
			kc := NewKsCalculator()

			// do calculation for each job.
			for job := range jobs {
				r1, r2 := job.r1, job.r2
				read1 := mated2Ref(r1)
				read2 := mated2Ref(r2)

				// double check if overlap.
				if r2.ReadLeft.Pos+len(read2) > r1.ReadLeft.Pos {
					// Determine overlap regions (in genome coordinate).
					start := maxInt(r1.ReadLeft.Pos, r2.ReadLeft.Pos)
					end := minInt(r1.ReadLeft.Pos+len(read1), r2.ReadLeft.Pos+len(read2))

					// Prepare profile and read sequences.
					profile := genome.PosProfile[start:end]
					nucl1 := read1[start-r1.ReadLeft.Pos : end-r1.ReadLeft.Pos]
					nucl2 := read2[start-r2.ReadLeft.Pos : end-r2.ReadLeft.Pos]

					// Subsitution profiling.
					subs := SubProfile(nucl1, nucl2, profile, pos)

					SubCorr(subs, cc, kc, maxl)
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

// Calculate correlations of substitutions in genomic sequences.
func CovGenomesGenomes(alignments []ncbiutils.SeqRecords, genome Genome, maxl, pos int) (kc *KsCalculator, cc *CovCalculator) {
	// Create job channel.
	jobs := make(chan ncbiutils.SeqRecords)
	go func() {
		defer close(jobs)
		for _, records := range alignments {
			jobs <- records
		}
	}()

	ncpu := runtime.GOMAXPROCS(0)

	// Create result channel.
	type result struct {
		cc *CovCalculator
		kc *KsCalculator
	}
	results := make(chan result)
	for i := 0; i < ncpu; i++ {
		go func() {
			cc := NewCovCalculator(maxl, true)
			kc := NewKsCalculator()
			acc := FindRefAcc(genome.Accession)
			for records := range jobs {
				refRecords := ncbiutils.SeqRecords{}
				for _, r := range records {
					acc2 := FindRefAcc(r.Genome)
					if acc == acc2 {
						refRecords = append(refRecords, r)
					}
				}

				for _, ref := range refRecords {
					// determine reference position profile.
					nucl := []byte{}
					for _, b := range ref.Nucl {
						if b != '-' {
							nucl = append(nucl, b)
						}
					}

					start := ref.Loc.From - 1
					end := ref.Loc.To - 3 // stop codon!
					var prof []byte
					if end > start {
						if end-start != len(nucl) {
							if end-start < len(nucl) {
								end = start + len(nucl)
							} else {
								Warn.Printf("Protein: %s profile length %d sequence length %d\n", ref.Id, end-start, len(nucl))
								continue
							}
						}
						prof = genome.PosProfile[start:end]
					} else {
						if len(genome.PosProfile)-start+end > len(nucl) {
							continue
						} else {
							end = len(nucl) + start - len(genome.PosProfile)
						}
						prof = genome.PosProfile[start:]
						prof = append(prof, genome.PosProfile[:end]...)
					}

					// compare substituations according to the profile.
					for j := 0; j < len(records); j++ {
						read1 := []byte{}
						for i, b := range ref.Nucl {
							if b != '-' {
								read1 = append(read1, records[j].Nucl[i])
							}
						}

						for k := j + 1; k < len(records); k++ {
							read2 := []byte{}
							for i, b := range ref.Nucl {
								if b != '-' {
									read2 = append(read2, records[k].Nucl[i])
								}
							}
							subs := SubProfile(read1, read2, prof, pos)
							SubCorr(subs, cc, kc, maxl)
						}
					}
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

// Calculate correlations of substituions in genomic sequences.
func CovGenomesGenome(alignments []ncbiutils.SeqRecords, genome Genome, maxl, pos int) (kc *KsCalculator, cc *CovCalculator) {
	// Create job channel.
	jobs := make(chan ncbiutils.SeqRecords)
	go func() {
		defer close(jobs)
		for _, records := range alignments {
			jobs <- records
		}
	}()

	ncpu := runtime.GOMAXPROCS(0)

	// Create result channel.
	type result struct {
		cc *CovCalculator
		kc *KsCalculator
	}
	results := make(chan result)
	for i := 0; i < ncpu; i++ {
		go func() {
			cc := NewCovCalculator(maxl, true)
			kc := NewKsCalculator()
			acc := FindRefAcc(genome.Accession)
			for records := range jobs {
				refRecords := ncbiutils.SeqRecords{}
				for _, r := range records {
					acc2 := FindRefAcc(r.Genome)
					if acc == acc2 {
						refRecords = append(refRecords, r)
					}
				}

				for _, ref := range refRecords {
					nucl := []byte{}
					for _, b := range ref.Nucl {
						if b != '-' {
							nucl = append(nucl, b)
						}
					}

					start := ref.Loc.From - 1
					end := ref.Loc.To - 3 // stop codon!
					var prof []byte
					if end > start {
						if end-start != len(nucl) {
							if end-start < len(nucl) {
								end = start + len(nucl)
							} else {
								Warn.Printf("Protein: %s profile length %d sequence length %d\n", ref.Id, end-start, len(nucl))
								continue
							}
						}
						prof = genome.PosProfile[start:end]
					} else {
						if len(genome.PosProfile)-start+end > len(nucl) {
							continue
						} else {
							end = len(nucl) + start - len(genome.PosProfile)
						}
						prof = genome.PosProfile[start:]
						prof = append(prof, genome.PosProfile[:end]...)
					}

					for _, r := range records {
						if ref.Genome != r.Genome {
							read := []byte{}
							for i, b := range ref.Nucl {
								if b != '-' {
									read = append(read, r.Nucl[i])
								}
							}

							subs := SubProfile(read, nucl, prof, pos)
							SubCorr(subs, cc, kc, maxl)
						}
					}
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
	subs := make([]float64, len(nucl))
	for i := 0; i < len(subs); i++ {
		p := profile[i]
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

// Calculate correlations of subsitutions.
func SubCorr(subs []float64, cc *CovCalculator, kc *KsCalculator, maxl int) {
	// Calculate kc and obtain position indices.
	posIndices := []int{}
	for k := 0; k < len(subs); k++ {
		if !math.IsNaN(subs[k]) {
			kc.Increment(subs[k])
			posIndices = append(posIndices, k)
			cc.Increment(0, subs[k], subs[k])
		}
	}
	// Calculate cc.
	for j := 0; j < len(posIndices); j++ {
		for k := j + 1; k < len(posIndices); k++ {
			l := posIndices[k] - posIndices[j]
			if l >= maxl {
				break
			} else {
				x, y := subs[posIndices[j]], subs[posIndices[k]]
				cc.Increment(l, x, y)
			}
		}
	}
}
