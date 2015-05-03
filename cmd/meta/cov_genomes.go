package main

import (
	"encoding/json"
	"flag"
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
	core      bool // whether to use core genomes.
	cmdConfig      // embed cmdConfig
}

func (cmd *cmdCovGenomes) Flags(fs *flag.FlagSet) *flag.FlagSet {
	fs = cmd.cmdConfig.Flags(fs)
	fs.BoolVar(&cmd.core, "core", false, "whether to use core genomes")
	return fs
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
		INFO.Printf("Total number of alignments: %d\n", len(alignments))
		var alns []seqrecord.SeqRecords
		if cmd.core {
			for _, aln := range alignments {
				m := make(map[string]bool)
				for _, rec := range aln {
					m[rec.Genome] = true
				}
				if len(m) == len(strains) {
					alns = append(alns, aln)
				}
			}
			INFO.Printf("Using core alignments: %d in total\n", len(alns))
		} else {
			alns = alignments
		}

		// If zero alignments, skip it.
		if len(alns) == 0 {
			WARN.Printf("%s has zero alignments.\n", prefix)
			continue
		}

		// For each position, do the calculation.
		for _, pos := range cmd.positions {
			cmd.RunOne(strains, alns, pos)
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

				covGenomesFuncs := []cov.GenomesOneFunc{
					cov.GenomesVsGenomeOne,
					cov.GenomesVsGenomesOne,
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

					if cmd.numBoot > 0 {
						results := cmd.covBoot(alignments, g, pos, covGenomesFunc)
						// Write result to files.
						filePrefix := fmt.Sprintf("%s_%s_pos%d", g.RefAcc(),
							funcType, pos)
						filePath := filepath.Join(*cmd.workspace, cmd.covOutBase, s.Path,
							filePrefix+"_boot.json")
						selected := []CovResult{}
						for _, res := range results {
							if !math.IsNaN(res.VarKs) {
								selected = append(selected, res)
							}
						}
						save2Jsons(selected, filePath)
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
func (cmd *cmdCovGenomes) Cov(records []seqrecord.SeqRecords, g genome.Genome, pos int, covFunc cov.GenomesOneFunc) (res CovResult) {
	kc, cc := cov.GenomesCalc(records, g, cmd.maxl, pos, covFunc)
	maxl := cmd.maxl
	res = createCovResult(kc, cc, maxl, pos)
	return
}

func (cmd *cmdCovGenomes) covBoot(records []seqrecord.SeqRecords, g genome.Genome, pos int, covFunc cov.GenomesOneFunc) (results []CovResult) {
	kcs, ccs := cov.GenomesBoot(records, g, cmd.maxl, pos, cmd.numBoot, covFunc)
	maxl := cmd.maxl
	for i := 0; i < len(kcs); i++ {
		kc, cc := kcs[i], ccs[i]
		res := createCovResult(kc, cc, maxl, pos)
		results = append(results, res)
	}
	return
}

func createCovResult(kc *cov.KsCalculator, cc *cov.CovCalculator, maxl, pos int) (res CovResult) {
	// Process and return a cov result.
	res.Ks = kc.Mean.GetResult()
	res.VarKs = kc.Var.GetResult()
	res.N = kc.Mean.GetN()

	// To use base distiance (step = 1) or codon distance (step = 3)
	var step, size int
	if pos == 0 {
		step = 1
		size = maxl
	} else {
		step = 3
		size = maxl / 3
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
