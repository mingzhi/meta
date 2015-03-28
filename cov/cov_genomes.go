package cov

import (
	"github.com/mingzhi/meta/genome"
	"github.com/mingzhi/ncbiftp/seqrecord"
	"log"
	"runtime"
)

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
			acc := g.RefAcc()
			for records := range jobs {
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
							continue
						} else {
							end = len(nucl) + start - len(g.PosProfile)
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
			acc := g.RefAcc()
			for records := range jobs {
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
							continue
						} else {
							end = len(nucl) + start - len(g.PosProfile)
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
