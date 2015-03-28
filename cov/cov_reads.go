package cov

import (
	"github.com/mingzhi/meta/genome"
	"github.com/mingzhi/meta/reads"
	"log"
	"runtime"
)

// Calculate correlation of substituions in reads,
// by comparing them to the reference genome.
// records: SamRecords;
// genome: Genome;
// maxl: max length of correlations;
// pos: positions to be calculated.
func ReadsVsGenome(matedReads reads.PairedEndReads, g genome.Genome, maxl, pos int) (kc *KsCalculator, cc *CovCalculator) {
	// Prepare jobs.
	type job struct {
		r reads.PairedEndRead
	}
	jobs := make(chan job)
	go func() {
		for _, r := range matedReads {
			sameRef := r.ReadLeft.Ref.Name() == r.ReadRight.Ref.Name()
			matchRef := genome.FindRefAcc(r.ReadLeft.Ref.Name()) == g.RefAcc()
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
				read := reads.MapMated2Ref(rec)

				if rec.ReadLeft.Pos+len(read) <= len(g.Seq) {
					start := rec.ReadLeft.Pos
					end := rec.ReadLeft.Pos + len(read)
					nucl := g.Seq[start:end]
					profile := g.PosProfile[start:end]
					subs := SubProfile(read, nucl, profile, pos)
					SubCorr(subs, cc, kc, maxl)
				} else {
					log.Printf("%d, %d, %d\n", rec.ReadLeft.Pos-1, rec.ReadLeft.Pos+len(read), len(g.PosProfile))
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
func ReadsVsReads(matedReads reads.PairedEndReads, g genome.Genome, maxl, pos int) (kc *KsCalculator, cc *CovCalculator) {
	ncpu := runtime.GOMAXPROCS(0)
	// Create job channel.
	type job struct {
		r1, r2 reads.PairedEndRead
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
				read1 := reads.MapMated2Ref(r1)
				read2 := reads.MapMated2Ref(r2)

				// double check if overlap.
				if r2.ReadLeft.Pos+len(read2) > r1.ReadLeft.Pos {
					// Determine overlap regions (in genome coordinate).
					start := maxInt(r1.ReadLeft.Pos, r2.ReadLeft.Pos)
					end := minInt(r1.ReadLeft.Pos+len(read1), r2.ReadLeft.Pos+len(read2))

					// Prepare profile and read sequences.
					profile := g.PosProfile[start:end]
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

// return max int
func maxInt(a, b int) int {
	if a > b {
		return a
	} else {
		return b
	}
}

// return min int
func minInt(a, b int) int {
	if a < b {
		return a
	} else {
		return b
	}
}
