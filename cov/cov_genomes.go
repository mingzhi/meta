package cov

import (
	"github.com/mingzhi/biogo/seq"
	"github.com/mingzhi/meta/genome"
	"github.com/mingzhi/ncbiftp/seqrecord"
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

// strip gaps.
func stripGaps(nucl []byte) []byte {
	return stripRefGaps(nucl, nucl)
}

// strip reference-gapped positions.
func stripRefGaps(ref []byte, read []byte) []byte {
	read1 := []byte{}
	for i, b := range ref {
		if isValidNucl(b) {
			read1 = append(read1, read[i])
		}
	}
	return read1
}

// get rightmost matched position profile.
func getProfile(rec seqrecord.SeqRecord, g genome.Genome) (prof []byte) {
	start := rec.Loc.From - 1
	end := rec.Loc.To

	var nucl []byte
	if end > start {
		nucl = g.Seq[start:end]
		prof = g.PosProfile[start:end]
		if rec.Loc.Strand == "-" {
			nucl = seq.Reverse(seq.Complement(nucl))
			prof = seq.Reverse(prof)
		}
	} else {
		nucl = g.Seq[start:]
		nucl = append(nucl, g.Seq[:end]...)
		prof = g.PosProfile[start:]
		prof = append(prof, g.PosProfile[:end]...)
	}

	read := stripGaps(rec.Nucl)
	l := 0
	for ; l < len(read) && l < len(nucl); l++ {
		if read[l] != nucl[l] {
			break
		}
	}

	prof = prof[:l]

	return
}

// calculate correlations.
func subMatrixCorr(subMatrix [][]float64, maxl int, c *Calculators) {
	for _, subs := range subMatrix {
		SubCorr(subs, c.TCov, c.Ks, maxl)
	}

	SubMatrixSCov(subMatrix, c.SCov, maxl)
	SubMatrixMCov(subMatrix, c.MCov, maxl)
	SubMatrixRCov(subMatrix, c.RCov, maxl)
}

func getRefRecords(records seqrecord.SeqRecords, g genome.Genome) (refRecords seqrecord.SeqRecords) {
	for _, r := range records {
		acc := genome.FindRefAcc(r.Genome)
		if acc == g.RefAcc() {
			refRecords = append(refRecords, r)
		}
	}
	return
}

type readPair struct {
	read1, read2 []byte
}

type readsVsFunc func(ref seqrecord.SeqRecord, records seqrecord.SeqRecords) (pairs []readPair)

func readsVsReads(ref seqrecord.SeqRecord, records seqrecord.SeqRecords) (pairs []readPair) {
	for i := 0; i < len(records); i++ {
		read1 := stripRefGaps(ref.Nucl, records[i].Nucl)
		for j := i + 1; j < len(records); j++ {
			read2 := stripRefGaps(ref.Nucl, records[j].Nucl)
			pairs = append(pairs, readPair{read1, read2})
		}
	}

	return
}

func readsVsRef(ref seqrecord.SeqRecord, records seqrecord.SeqRecords) (pairs []readPair) {
	read1 := stripGaps(ref.Nucl)
	for i := 0; i < len(records); i++ {
		if records[i].Genome != ref.Genome {
			read2 := stripRefGaps(ref.Nucl, records[i].Nucl)
			pairs = append(pairs, readPair{read1, read2})
		}
	}

	return
}

func genomeOne(records seqrecord.SeqRecords, g genome.Genome, maxl, pos int, c *Calculators, f readsVsFunc) {
	refRecords := getRefRecords(records, g)
	for _, ref := range refRecords {
		// determine reference position profile.
		prof := getProfile(ref, g)
		if len(prof) < len(ref.Nucl)/2 {
			log.Printf("Discard sequence record:\n")
			log.Printf("%s, %s, %v\n", ref.Id, ref.Genome, ref.Loc)
			continue
		}

		pairs := f(ref, records)
		subMatrix := [][]float64{}
		for _, pair := range pairs {
			read1 := pair.read1[:len(prof)]
			read2 := pair.read2[:len(prof)]
			subs := SubProfile(read1, read2, prof, pos)
			subMatrix = append(subMatrix, subs)
		}
		subMatrixCorr(subMatrix, maxl, c)
	}
}

func GenomesVsGenomesOne(records seqrecord.SeqRecords, g genome.Genome, maxl, pos int, c *Calculators) {
	f := readsVsReads
	genomeOne(records, g, maxl, pos, c, f)
}

func GenomesVsGenomeOne(records seqrecord.SeqRecords, g genome.Genome, maxl, pos int, c *Calculators) {
	f := readsVsRef
	genomeOne(records, g, maxl, pos, c, f)
}
