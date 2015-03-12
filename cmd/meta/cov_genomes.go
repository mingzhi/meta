package main

import (
	"encoding/json"
	"fmt"
	"github.com/mingzhi/meta"
	"github.com/mingzhi/ncbiutils"
	"math"
	"os"
	"path/filepath"
)

// Command to calculate correlations of substituions
// in reference genomes.
type cmdCovGenomes struct {
	cmdConfig // embed cmdConfig
}

// Run command.
func (cmd *cmdCovGenomes) Run(args []string) {
	// Parse config and settings.
	cmd.ParseConfig()
	MakeDir(cmd.covOutBase)

	// Read strain information.
	strainFilePath := filepath.Join(*cmd.workspace, cmd.strainFileName)
	strains := meta.ReadStrains(strainFilePath)

	// Make position profiles.
	meta.GenomePosProfiling(strains, cmd.refBase)

	// Read alignments.
	alignments := cmd.ReadAlignments()

	// For each strain.
	for _, s := range strains {
		for _, genome := range s.Genomes {
			if isChromosome(genome.Replicon) {
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
					res := cmd.Cov(alignments, genome, pos)
					// Write result to files.
					filePrefix := fmt.Sprintf("%s_%s_pos%d", s.Path,
						"Cov_Genomes_vs_Genomes", pos)
					filePath := filepath.Join(*cmd.workspace, cmd.covOutBase,
						filePrefix+".json")
					save2Json(res, filePath)
				}
			}
		}
	}
}

// Load alignments.
func (cmd *cmdCovGenomes) ReadAlignments() (alns []ncbiutils.SeqRecords) {
	fileName := cmd.prefix + "_orthologs_aligned.json"
	filePath := filepath.Join(*cmd.workspace, cmd.orthoOutBase,
		fileName)
	r, err := os.Open(filePath)
	if err != nil {
		ERROR.Fatalln(err)
	}
	defer r.Close()

	decoder := json.NewDecoder(r)
	err = decoder.Decode(&alns)
	if err != nil {
		ERROR.Fatalln(err)
	}

	return
}

// Calculate covariances for records.
func (cmd *cmdCovGenomes) Cov(records []ncbiutils.SeqRecords,
	genome meta.Genome, pos int) (res CovResult) {
	kc, cc := meta.CovGenomes(records, genome, cmd.maxl, pos)

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
