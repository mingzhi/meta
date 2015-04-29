package cov

import (
	"github.com/mingzhi/meta/genome"
	"github.com/mingzhi/ncbiftp/seqrecord"
	"log"
	"math/rand"
	"runtime"
)

type GenomesOneFunc func(records seqrecord.SeqRecords, g genome.Genome, maxl, pos int, kc *KsCalculator, cc *CovCalculator)

func GenomesBoot(alignments []seqrecord.SeqRecords, g genome.Genome, maxl, pos, numBoot int, oneFunc GenomesOneFunc) ([]*KsCalculator, []*CovCalculator) {
	biasCorrection := true
	// Create job channel.
	jobs := make(chan seqrecord.SeqRecords)
	go func() {
		defer close(jobs)
		for _, records := range alignments {
			jobs <- records
		}
	}()

	ncpu := runtime.GOMAXPROCS(0)

	type result struct {
		cc *CovCalculator
		kc *KsCalculator
	}

	done := make(chan bool)

	resultChan := make(chan result)
	for i := 0; i < ncpu; i++ {
		go func() {
			for records := range jobs {
				kc := NewKsCalculator()
				cc := NewCovCalculator(maxl, biasCorrection)
				oneFunc(records, g, maxl, pos, kc, cc)
				res := result{cc: cc, kc: kc}
				resultChan <- res
			}
			done <- true
		}()
	}

	go func() {
		defer close(resultChan)
		for i := 0; i < ncpu; i++ {
			<-done
		}
	}()

	results := []result{}
	for res := range resultChan {
		results = append(results, res)
	}

	// bootstrapping
	bootJobs := make(chan []int)
	go func() {
		defer close(bootJobs)
		for i := 0; i < numBoot; i++ {
			sample := make([]int, len(results))
			for j := 0; j < len(sample); j++ {
				sample[j] = rand.Intn(len(results))
			}
			bootJobs <- sample
		}
	}()

	bootResultChan := make(chan result)
	for i := 0; i < ncpu; i++ {
		go func() {
			for sample := range bootJobs {
				kc := NewKsCalculator()
				cc := NewCovCalculator(maxl, biasCorrection)
				for j := 0; j < len(sample); j++ {
					index := sample[j]
					res := results[index]
					kc.Append(res.kc)
					cc.Append(res.cc)
				}
				bootResultChan <- result{cc: cc, kc: kc}
			}
			done <- true
		}()
	}

	go func() {
		defer close(bootResultChan)
		for i := 0; i < ncpu; i++ {
			<-done
		}
	}()

	kcs := []*KsCalculator{}
	ccs := []*CovCalculator{}
	for res := range bootResultChan {
		kcs = append(kcs, res.kc)
		ccs = append(ccs, res.cc)
	}

	return kcs, ccs
}

func GenomesCalc(alignments []seqrecord.SeqRecords, g genome.Genome, maxl, pos int, covFunc GenomesOneFunc) (kc *KsCalculator, cc *CovCalculator) {
	jobs := make(chan seqrecord.SeqRecords)
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
			for records := range jobs {
				covFunc(records, g, maxl, pos, kc, cc)
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

type GenomesFunc func(alignments []seqrecord.SeqRecords, g genome.Genome, maxl, pos int) (kc *KsCalculator, cc *CovCalculator)

// Calculate correlations of substitutions in genomic sequences.
func GenomesVsGenomes(alignments []seqrecord.SeqRecords, g genome.Genome, maxl, pos int) (kc *KsCalculator, cc *CovCalculator) {
	// Create job channel.
	jobs := make(chan seqrecord.SeqRecords)
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
			for records := range jobs {
				GenomesVsGenomesOne(records, g, maxl, pos, kc, cc)
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

func GenomesVsGenomesOne(records seqrecord.SeqRecords, g genome.Genome, maxl, pos int,
	kc *KsCalculator, cc *CovCalculator) {
	acc := g.RefAcc()
	refRecords := seqrecord.SeqRecords{}
	for _, r := range records {
		acc2 := genome.FindRefAcc(r.Genome)
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
					log.Printf("Protein: %s profile length %d sequence length %d\n", ref.Id, end-start, len(nucl))
					continue
				}
			}
			prof = g.PosProfile[start:end]
		} else {
			if len(g.PosProfile)-start+end > len(nucl) {
				// [TODO]
				continue
			} else {
				end = len(nucl) + start - len(g.PosProfile)
				// [TODO]
				if end <= 0 {
					continue
				}
			}
			prof = g.PosProfile[start:]
			prof = append(prof, g.PosProfile[:end]...)
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

// Calculate correlations of substituions in genomic sequences.
func GenomesVsGenome(alignments []seqrecord.SeqRecords, g genome.Genome, maxl, pos int) (kc *KsCalculator, cc *CovCalculator) {
	// Create job channel.
	jobs := make(chan seqrecord.SeqRecords)
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
			for records := range jobs {
				GenomesVsGenomeOne(records, g, maxl, pos, kc, cc)
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

func GenomesVsGenomeOne(records seqrecord.SeqRecords, g genome.Genome, maxl, pos int,
	kc *KsCalculator, cc *CovCalculator) {
	acc := g.RefAcc()
	refRecords := seqrecord.SeqRecords{}
	for _, r := range records {
		acc2 := genome.FindRefAcc(r.Genome)
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
					log.Printf("Protein: %s profile length %d sequence length %d\n", ref.Id, end-start, len(nucl))
					continue
				}
			}
			prof = g.PosProfile[start:end]
		} else {
			if len(g.PosProfile)-start+end > len(nucl) {
				// [TODO]
				continue
			} else {
				end = len(nucl) + start - len(g.PosProfile)
				if end <= 0 {
					// [TODO]
					continue
				}
			}
			prof = g.PosProfile[start:]
			prof = append(prof, g.PosProfile[:end]...)
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
