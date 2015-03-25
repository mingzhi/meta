package main

import (
	"fmt"
	"github.com/mingzhi/meta"
	"math"
	"os"
	"path/filepath"
	"sort"
	"strings"
)

type covReadsFunc func(records meta.PairedEndReads,
	genome meta.Genome, maxl, pos int) (kc *meta.KsCalculator, cc *meta.CovCalculator)

// Command to calculate correlations for mapped reads to reference genomes.
type cmdCovReads struct {
	cmdConfig // embedded cmdConfig.

	covFunc covReadsFunc // cov calculate function.
}

// Run command.
func (cmd *cmdCovReads) Run(args []string) {
	// Parse config and settings.
	cmd.ParseConfig()
	// Load species:strains map.
	cmd.LoadSpeciesMap()
	// Make cov output diretory.
	MakeDir(filepath.Join(*cmd.workspace, cmd.covOutBase))

	for _, strains := range cmd.speciesMap {
		// For each strain.
		for _, s := range strains {
			// Make specific output folder.
			MakeDir(filepath.Join(*cmd.workspace, cmd.covOutBase, s.Path))
			for _, genome := range s.Genomes {
				// [TODO] Now we only work on chromosome genomes.
				if isChromosome(genome.Replicon) {
					// Clear RefSeq accession (remove .version and such).
					acc := meta.FindRefAcc(genome.Accession)
					// Read records of reads from a "sam" file.
					samFileName := acc + bowtiedSamAppendix
					samFilePath := filepath.Join(*cmd.workspace, cmd.samOutBase, s.Path, samFileName)
					// Check if the "sam" file exists.
					if isSamFileExist(samFilePath) {
						_, records := meta.ReadSamFile(samFilePath)

						if len(records) == 0 {
							WARN.Printf("%s,%s has zero records\n", s.Path, genome.Accession)
						} else {
							// Read position profile for the genome.
							posFileName := meta.FindRefAcc(genome.Accession) + ".pos"
							posFilePath := filepath.Join(cmd.refBase, s.Path, posFileName)
							genome.PosProfile = meta.ReadPosProfile(posFilePath)
							// Read sequence for the genome.
							fnaFileName := meta.FindRefAcc(genome.Accession) + ".fna"
							fnaFilePath := filepath.Join(cmd.refBase, s.Path, fnaFileName)
							genome.Seq = meta.ReadFasta(fnaFilePath).Seq

							// paired end reads and sorted.
							matedReads := meta.GetPairedEndReads(records)
							INFO.Printf("%s has %d mated reads.\n", samFilePath, len(matedReads))

							for _, funcName := range cmd.covReadsFunctions {
								// Assign cov read function.
								switch funcName {
								case "Cov_Reads_vs_Reads":
									cmd.covFunc = meta.CovReadsReads
									sort.Sort(meta.ByRightCoordinatePairedEndReads{matedReads})
								case "Cov_Reads_vs_Genome":
									cmd.covFunc = meta.CovReadsGenome
								default:
									continue
								}

								// Calculate correlations at each position.
								for _, pos := range cmd.positions {
									res := cmd.Cov(matedReads, genome, pos)
									// Write result to files.
									filePrefix := fmt.Sprintf("%s_%s_pos%d", acc,
										funcName, pos)
									filePath := filepath.Join(*cmd.workspace, cmd.covOutBase, s.Path,
										filePrefix+".json")
									if !math.IsNaN(res.VarKs) {
										res.NReads = len(matedReads)
										save2Json(res, filePath)
									} else {
										WARN.Printf("%s: VarKs: NaN\n", filePath)
									}

								}
							}

						}
					} else {
						WARN.Printf("Cannot find sam file: %s\n", samFilePath)
					}
				}
			}
		}
	}

}

// Calculate covariance for records.
func (cmd *cmdCovReads) Cov(records meta.PairedEndReads,
	genome meta.Genome, pos int) (res CovResult) {

	kc, cc := cmd.covFunc(records, genome, cmd.maxl, pos)

	// Process and return a cov result.
	res.Ks = kc.Mean.GetResult()
	res.VarKs = kc.Var.GetResult()
	res.N = kc.Mean.GetN()

	// To use base distiance (step = 1) or codon distance (step = 3)
	var step, size int
	if pos == 0 {
		step = 1
		size = cmd.maxl
	} else {
		step = 3
		size = cmd.maxl / 3
	}

	for i := 0; i < size; i++ {
		index := step * i
		v := cc.GetResult(index)
		n := cc.GetN(index)
		if !math.IsNaN(v) {
			res.CtIndices = append(res.CtIndices, i)
			res.Ct = append(res.Ct, v)
			res.CtN = append(res.CtN, n)
		}
	}

	return
}

// Check if it is a chromosome,
// by simply searching "chromosome" keyword.
func isChromosome(replicon string) (is bool) {
	if strings.Contains(replicon, "chromosome") {
		is = true
	} else {
		is = false
	}
	return
}

// Check if the sam file exists.
func isSamFileExist(filePath string) (isExist bool) {
	_, err := os.Stat(filePath)
	if err != nil {
		if os.IsNotExist(err) {
			WARN.Printf("%s does not exist!\n", filePath)
			isExist = false
		} else {
			ERROR.Fatalln(err)
		}
	} else {
		isExist = true
	}

	return
}
