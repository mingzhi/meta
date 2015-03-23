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
	cmd.LoadSpeciesMap()
	MakeDir(filepath.Join(*cmd.workspace, cmd.covOutBase))

	for prefix, strains := range cmd.speciesMap {
		// Make position profiles.
		meta.GenomePosProfiling(strains, cmd.refBase)

		// Read alignments.
		alignments := cmd.ReadAlignments(prefix)

		for _, pos := range cmd.positions {
			cmd.RunOne(strains, alignments, pos)
		}
	}

}

func (cmd *cmdCovGenomes) RunOne(strains []meta.Strain, alignments []ncbiutils.SeqRecords, pos int) {
	// For each strain (genome), create a job.
	type job struct {
		strain meta.Strain
		genome meta.Genome
	}
	jobs := make(chan job)
	go func() {
		defer close(jobs)
		for _, strain := range strains {
			MakeDir(filepath.Join(*cmd.workspace, cmd.covOutBase, strain.Path))
			for _, genome := range strain.Genomes {
				if isChromosome(genome.Replicon) {
					jobs <- job{strain, genome}
				}
			}
		}
	}()

	done := make(chan bool)
	for i := 0; i < *cmd.ncpu; i++ {
		go func() {
			for job := range jobs {
				strain := job.strain
				genome := job.genome
				// Read position profile for the genome.
				posFileName := meta.FindRefAcc(genome.Accession) + ".pos"
				posFilePath := filepath.Join(cmd.refBase, strain.Path, posFileName)
				genome.PosProfile = meta.ReadPosProfile(posFilePath)
				// Read sequence for the genome.
				fnaFileName := meta.FindRefAcc(genome.Accession) + ".fna"
				fnaFilePath := filepath.Join(cmd.refBase, strain.Path, fnaFileName)
				genome.Seq = meta.ReadFasta(fnaFilePath).Seq

				covGenomesFuncs := []meta.CovGenomesFunc{
					meta.CovGenomesGenome,
					meta.CovGenomesGenomes,
				}

				covGenomesFuncNames := []string{
					"Cov_Genomes_vs_Genome",
					"Cov_Genomes_vs_Genomes",
				}

				acc := meta.FindRefAcc(genome.Accession)

				for j, covGenomesFunc := range covGenomesFuncs {
					funcType := covGenomesFuncNames[j]
					res := cmd.Cov(alignments, genome, pos, covGenomesFunc)
					// Write result to files.
					filePrefix := fmt.Sprintf("%s_%s_pos%d", acc,
						funcType, pos)
					filePath := filepath.Join(*cmd.workspace, cmd.covOutBase, strain.Path,
						filePrefix+".json")
					save2Json(res, filePath)
				}

			}
			done <- true
		}()
	}

	// Waiting for the job done.
	for i := 0; i < *cmd.ncpu; i++ {
		<-done
	}
}

// Load alignments.
func (cmd *cmdCovGenomes) ReadAlignments(prefix string) (alns []ncbiutils.SeqRecords) {
	fileName := prefix + "_orthologs_aligned.json"
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
	genome meta.Genome, pos int, covGenomesFunc meta.CovGenomesFunc) (res CovResult) {
	kc, cc := covGenomesFunc(records, genome, cmd.maxl, pos)

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
