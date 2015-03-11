package meta

import (
	"github.com/mingzhi/biogo/seq"
	"github.com/mingzhi/ncbiutils"
	"log"
	"os"
	"path/filepath"
	"runtime"
	"strings"
)

const (
	FirstPos byte = 1 << iota
	SecondPos
	ThirdPos
	FourFold
)

// Generate a position profile for a genome.
// strain: Strain contains the path and genome name.
// dir: folder containing genome sequences.
// TODO: profile boundary regions.
func GenomePosProfiling(strains []Strain, dir string) {
	// Prepare genetic code table map.
	gcMap := ncbiutils.GeneticCodes()

	// Prepare job queue.
	ncpu := runtime.GOMAXPROCS(0)
	jobs := make(chan Strain)
	go func() {
		for _, s := range strains {
			jobs <- s
		}
		close(jobs)
	}()

	done := make(chan bool)
	for i := 0; i < ncpu; i++ {
		go func() {
			for s := range jobs {
				sDir := filepath.Join(dir, s.Path) // strain folder
				// Genetic code table, we need it for 4-fold codons.
				gc := gcMap[s.GeneticCode]
				for _, g := range s.Genomes {
					if !strings.Contains(g.Replicon, "chromosome") {
						continue
					}
					// Read genome sequence from fna file.
					acc := FindRefAcc(g.Accession)
					fnaFileName := acc + ".fna"
					fnaFilePath := filepath.Join(sDir, fnaFileName)
					genomeSeq := readFasta(fnaFilePath)[0].Seq

					// Initalize position profile.
					profile := make([]byte, len(genomeSeq))

					// Read protein features from ptt file.
					pttFileName := acc + ".ptt"
					pttFilePath := filepath.Join(sDir, pttFileName)
					pttFile := ncbiutils.NewPttFile(pttFilePath)
					ptts := pttFile.ReadAll()

					// For each gene (a ptt), record codon positions.
					for _, ptt := range ptts {
						// Prepare nucleotide sequence,
						// we need it for determine 4-fold codons.
						var nucl []byte
						if ptt.Loc.To >= ptt.Loc.From {
							nucl = genomeSeq[ptt.Loc.From-1 : ptt.Loc.To]
						} else {
							// skip genes across boundary.
							continue
						}

						if ptt.Loc.Strand == "-" {
							nucl = seq.Complement(seq.Reverse(nucl))
						}

						prof := make([]byte, len(nucl))
						for j, _ := range nucl {
							switch (j + 1) % 3 {
							case 1:
								prof[j] = FirstPos
							case 2:
								prof[j] = SecondPos
							case 0:
								codon := nucl[j-2 : j+1]
								if gc.FFCodons[string(codon)] {
									prof[j] = FourFold
								} else {
									prof[j] = ThirdPos
								}
							}
						}

						if ptt.Loc.Strand == "-" {
							prof = seq.Reverse(prof)
						}

						for j, p := range prof {
							profile[ptt.Loc.From-1+j] = p
						}
					}

					// Save genome profile to file.
					fileName := strings.Replace(fnaFilePath, "fna", "pos", -1)
					writePosProfile(fileName, profile)
				}
			}
			done <- true
		}()
	}

	for i := 0; i < ncpu; i++ {
		<-done
	}
}

func readFasta(fileName string) []*seq.Sequence {
	f, err := os.Open(fileName)
	if err != nil {
		log.Panic(err)
	}
	defer f.Close()

	fastaRd := seq.NewFastaReader(f)
	seqs, err := fastaRd.ReadAll()
	if err != nil {
		log.Panic(err)
	}

	return seqs
}

func writePosProfile(fileName string, profile []byte) {
	w, err := os.Create(fileName)
	if err != nil {
		log.Panic(err)
	}
	defer w.Close()

	w.Write(profile)
}
