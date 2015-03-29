package main

import (
	"encoding/json"
	"fmt"
	"github.com/mingzhi/meta/cov"
	"github.com/mingzhi/meta/genome"
	"github.com/mingzhi/meta/strain"
	"github.com/mingzhi/ncbiftp/seqrecord"
	"math"
	"os"
	"path/filepath"
)

// Command to calculate correlations of substituions
// in reference genomes.
type cmdCovGenomes struct {
	cmdConfig // embed cmdConfig
}

func (cmd *cmdCovGenomes) Init() {
	// Parse config and settings.
	cmd.ParseConfig()
	// Load species map.
	cmd.LoadSpeciesMap()
	// Make output directory.
	MakeDir(filepath.Join(*cmd.workspace, cmd.covOutBase))
	// Check profile positions.
	if len(cmd.positions) == 0 {
		WARN.Println("Use default position: 4!")
		cmd.positions = append(cmd.positions, 4)
	}
}

// Run command.
func (cmd *cmdCovGenomes) Run(args []string) {
	// Initialize.
	cmd.Init()

	// For each species in the species map.
	for prefix, strains := range cmd.speciesMap {
		// Read alignments.
		alignments := cmd.ReadAlignments(prefix)

		// If zero alignments, skip it.
		if len(alignments) == 0 {
			WARN.Printf("%s has zero alignments.\n", prefix)
			continue
		}

		// For each position, do the calculation.
		for _, pos := range cmd.positions {
			cmd.RunOne(strains, alignments, pos)
		}
	}

}

func (cmd *cmdCovGenomes) RunOne(strains []strain.Strain, alignments []seqrecord.SeqRecords, pos int) {
	// For each strain (genome), create a job.
	type job struct {
		strain strain.Strain
		genome genome.Genome
	}
	jobs := make(chan job)
	go func() {
		defer close(jobs)
		for _, strain := range strains {
			// Create output directory.
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
				s := job.strain
				g := job.genome
				// base folder of the strain.
				base := filepath.Join(cmd.refBase, s.Path)
				genome.LoadFna(&g, base)
				genome.LoadProfile(&g, base)

				covGenomesFuncs := []cov.GenomesFunc{
					cov.GenomesVsGenome,
					cov.GenomesVsGenomes,
				}

				covGenomesFuncNames := []string{
					"Cov_Genomes_vs_Genome",
					"Cov_Genomes_vs_Genomes",
				}

				for j, covGenomesFunc := range covGenomesFuncs {
					funcType := covGenomesFuncNames[j]
					res := cmd.Cov(alignments, g, pos, covGenomesFunc)
					// Write result to files.
					filePrefix := fmt.Sprintf("%s_%s_pos%d", g.RefAcc(),
						funcType, pos)
					filePath := filepath.Join(*cmd.workspace, cmd.covOutBase, s.Path,
						filePrefix+".json")
					if !math.IsNaN(res.VarKs) {
						save2Json(res, filePath)
					} else {
						WARN.Printf("%s: VarKs: NaN\n", filePath)
					}
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
func (cmd *cmdCovGenomes) ReadAlignments(prefix string) (alns []seqrecord.SeqRecords) {
	fileName := prefix + "_orthologs_aligned.json"
	filePath := filepath.Join(*cmd.workspace, cmd.orthoOutBase,
		fileName)
	r, err := os.Open(filePath)
	if err != nil {
		WARN.Println(err)
		return
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
func (cmd *cmdCovGenomes) Cov(records []seqrecord.SeqRecords,
	g genome.Genome, pos int, covGenomesFunc cov.GenomesFunc) (res CovResult) {
	kc, cc := covGenomesFunc(records, g, cmd.maxl, pos)

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
