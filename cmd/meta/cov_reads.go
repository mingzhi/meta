package main

import (
	"fmt"
	"github.com/mingzhi/meta/cov"
	"github.com/mingzhi/meta/genome"
	"github.com/mingzhi/meta/reads"
	"math"
	"os"
	"path/filepath"
	"sort"
	"strings"
)

type covReadsFunc func(records reads.PairedEndReads,
	g genome.Genome, maxl, pos int) (kc *cov.KsCalculator, cc *cov.CovCalculator)

// Command to calculate correlations for mapped reads to reference genomes.
type cmdCovReads struct {
	cmdConfig // embedded cmdConfig.

	covFunc covReadsFunc // cov calculate function.
}

func (cmd *cmdCovReads) Init() {
	if len(cmd.positions) == 0 {
		WARN.Println("Use default position: 4!")
	}
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
			for _, g := range s.Genomes {
				// [TODO] Now we only work on chromosome genomes.
				if isChromosome(g.Replicon) {
					// Read records of reads from a "sam" file.
					samFileName := g.RefAcc() + bowtiedSamAppendix
					samFilePath := filepath.Join(*cmd.workspace, cmd.samOutBase, s.Path, samFileName)
					// Check if the "sam" file exists.
					if isSamFileExist(samFilePath) {
						_, records := reads.ReadSamFile(samFilePath)

						if len(records) == 0 {
							WARN.Printf("%s,%s has zero records\n", s.Path, g.RefAcc())
						} else {
							// Read position profile for the genome.
							// base folder of the strain.
							base := filepath.Join(cmd.refBase, s.Path)
							genome.LoadFna(&g, base)
							genome.LoadProfile(&g, base)

							// paired end reads and sorted.
							matedReads := reads.GetPairedEndReads(records)
							for _, funcName := range cmd.covReadsFuncs {
								// Assign cov read function.
								switch funcName {
								case "Cov_Reads_vs_Reads":
									cmd.covFunc = cov.ReadsVsReads
									sort.Sort(reads.ByRightCoordinatePairedEndReads{matedReads})
								case "Cov_Reads_vs_Genome":
									cmd.covFunc = cov.ReadsVsGenome
								default:
									continue
								}

								// Calculate correlations at each position.
								for _, pos := range cmd.positions {
									res := cmd.Cov(matedReads, g, pos)
									// Write result to files.
									filePrefix := fmt.Sprintf("%s_%s_pos%d", g.RefAcc(),
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
func (cmd *cmdCovReads) Cov(records reads.PairedEndReads,
	g genome.Genome, pos int) (res CovResult) {

	kc, cc := cmd.covFunc(records, g, cmd.maxl, pos)

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
