package cov

import (
	"github.com/mingzhi/meta/genome"
	"github.com/mingzhi/ncbiftp/seqrecord"
	"github.com/mingzhi/ncbiftp/taxonomy"
	"log"
	"math/rand"
	"runtime"
)

type GenomesOneFunc func(records seqrecord.SeqRecords, g genome.Genome, maxl, pos int, c *Calculators)

func GenomesCalc(alignments []seqrecord.SeqRecords, g genome.Genome, maxl, pos int, oneFunc GenomesOneFunc) []*Calculators {
	return genomeCalc(alignments, g, maxl, pos, oneFunc)
}

func GenomesBoot(alignments []seqrecord.SeqRecords, g genome.Genome, maxl, pos, numBoot int, oneFunc GenomesOneFunc) <-chan []*Calculators {
	results := genomeCalc(alignments, g, maxl, pos, oneFunc)

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

	ncpu := runtime.GOMAXPROCS(0)
	bootResultChan := make(chan []*Calculators)
	done := make(chan bool)
	for i := 0; i < ncpu; i++ {
		go func() {
			for sample := range bootJobs {
				bootResults := []*Calculators{}
				for j := 0; j < len(sample); j++ {
					index := sample[j]
					res := results[index]
					bootResults = append(bootResults, res)
				}
				bootResultChan <- bootResults
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

	return bootResultChan
}

func genomeCalc(alignments []seqrecord.SeqRecords, g genome.Genome, maxl, pos int, oneFunc GenomesOneFunc) []*Calculators {
	biasCorrection := false
	// Create job channel.
	jobs := make(chan seqrecord.SeqRecords)
	go func() {
		defer close(jobs)
		for _, records := range alignments {
			jobs <- records
		}
	}()

	ncpu := runtime.GOMAXPROCS(0)

	done := make(chan bool)

	resultChan := make(chan *Calculators)
	for i := 0; i < ncpu; i++ {
		go func() {
			for records := range jobs {
				c := NewCalculators(maxl, biasCorrection)
				oneFunc(records, g, maxl, pos, c)
				resultChan <- c
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

	results := []*Calculators{}
	for res := range resultChan {
		results = append(results, res)
	}

	return results
}

func fourFoldProfile(records seqrecord.SeqRecords, gc *taxonomy.GeneticCode) []bool {
	n := len(records[0].Nucl) / 3
	profile := make([]bool, n)
	ffcondon := gc.FFCodons
	for k := 0; k < n; k++ {
		profile[k] = true
		codonA := string(records[0].Nucl[3*k : 3*k+3])
		aa := gc.Table[codonA]
		for i := 0; i < len(records); i++ {
			codon := string(records[i].Nucl[3*k : 3*k+3])
			ab := gc.Table[codon]
			if aa == '-' || aa != ab || !ffcondon[codon] {
				profile[k] = false
				break
			}
		}
	}

	profile1 := make([]bool, len(records[0].Nucl))
	for i := 0; i < n; i++ {
		profile1[3*i] = profile[i]
	}

	return profile1
}

var GCMAP map[string]*taxonomy.GeneticCode

func GenomesVsGenomesOne(records seqrecord.SeqRecords, g genome.Genome, maxl, pos int, c *Calculators) {
	acc := g.RefAcc()
	refRecords := seqrecord.SeqRecords{}
	for _, r := range records {
		acc2 := genome.FindRefAcc(r.Genome)
		if acc == acc2 {
			refRecords = append(refRecords, r)
		}
	}

	if GCMAP == nil {
		GCMAP = taxonomy.GeneticCodes()
	}

	gc := GCMAP[records[0].Code]

	ffProfile := fourFoldProfile(records, gc)

	for _, ref := range refRecords {
		// determine reference position profile.
		nucl := []byte{}
		ff := []bool{}
		for i, b := range ref.Nucl {
			if b != '-' {
				nucl = append(nucl, b)
				ff = append(ff, ffProfile[i])
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

		for i := 0; i < len(prof); i++ {
			if prof[i] == genome.FourFold {
				if !ff[i] {
					prof[i] = genome.ThirdPos
				}
			}
		}

		// compare substituations according to the profile.
		subMatrix := [][]float64{}
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
				SubCorr(subs, c.TCov, c.Ks, maxl)
				subMatrix = append(subMatrix, subs)
			}
		}
		SubMatrixSCov(subMatrix, c.SCov, maxl)
		SubMatrixMCov(subMatrix, c.MCov, maxl)
		SubMatrixRCov(subMatrix, c.RCov, maxl)
	}
}

func GenomesVsGenomeOne(records seqrecord.SeqRecords, g genome.Genome, maxl, pos int, c *Calculators) {
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

		subMatrix := [][]float64{}
		for _, r := range records {
			if ref.Genome != r.Genome {
				read := []byte{}
				for i, b := range ref.Nucl {
					if b != '-' {
						read = append(read, r.Nucl[i])
					}
				}

				subs := SubProfile(read, nucl, prof, pos)
				SubCorr(subs, c.TCov, c.Ks, maxl)
				subMatrix = append(subMatrix, subs)
			}
		}

		SubMatrixSCov(subMatrix, c.SCov, maxl)
		SubMatrixMCov(subMatrix, c.MCov, maxl)
		SubMatrixRCov(subMatrix, c.RCov, maxl)
	}

}
