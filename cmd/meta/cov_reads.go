package main

import (
	"encoding/json"
	"fmt"
	"github.com/mingzhi/meta"
	"math"
	"os"
	"path/filepath"
	"strings"
)

type CovResult struct {
	Ks, VarKs float64
	Ct        []float64
	CtIndices []int
	CtN       []int
	N         int
}

type covReadsFunc func(records meta.SamRecords,
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
	MakeDir(cmd.covOutBase)

	// Assign cov read function.
	switch cmd.covReadsFuncName {
	case "CovReads":
		cmd.covFunc = meta.CovReads
	default:
		cmd.covFunc = meta.CovGenome
	}

	// Read strain information.
	strainFilePath := filepath.Join(*cmd.workspace, cmd.strainFileName)
	strains := meta.ReadStrains(strainFilePath)

	// Make position profiles.
	meta.GenomePosProfiling(strains, cmd.refBase)

	// For each strain.
	for _, s := range strains {
		for _, genome := range s.Genomes {
			// We only work on chromosome genomes.
			if isChromosome(genome.Replicon) {
				// Read records of reads from a "sam" file.
				samFileName := s.Path + ".sam"
				samFilePath := filepath.Join(cmd.samOutBase, samFileName)
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

						// Calculate correlations at each position.
						for _, pos := range cmd.positions {
							res := cmd.Cov(records, genome, pos)
							// Write result to files.
							filePrefix := fmt.Sprintf("%s_%s_pos%d", s.Path,
								cmd.covReadsFuncName, pos)
							filePath := filepath.Join(*cmd.workspace, cmd.covOutBase,
								filePrefix+".json")
							save2Json(res, filePath)

						}
					}
				}
			}
		}
	}
}

// Calculate covariance for records.
func (cmd *cmdCovReads) Cov(records meta.SamRecords,
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

// Save cov result to a json file.
func save2Json(cr CovResult, fileName string) {
	f, err := os.Create(fileName)
	if err != nil {
		ERROR.Fatalln(err)
	}
	defer f.Close()

	ec := json.NewEncoder(f)
	if err := ec.Encode(cr); err != nil {
		ERROR.Fatalln(err)
	}
}
