package main

import (
	"compress/flate"
	"encoding/json"
	"flag"
	"fmt"
	"github.com/mingzhi/meta/cov"
	"github.com/mingzhi/meta/genome"
	"github.com/mingzhi/meta/strain"
	"github.com/mingzhi/ncbiftp/seqrecord"
	"log"
	"math"
	"math/rand"
	"os"
	"path/filepath"
	"strings"
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

	type job struct {
		strains    []strain.Strain
		alignments []seqrecord.SeqRecords
		position   int
		name       string
	}

	jobChan := make(chan job)
	go func() {
		defer close(jobChan)
		for prefix, strains := range cmd.speciesMap {
			for _, pos := range cmd.positions {
				p := prefix
				if pos == 0 {
					appendix := "expanded"
					p = strings.Join([]string{prefix, appendix}, "_")
				}
				alignments := cmd.ReadAlignments(p)
				cores, disps := separeteAlignments(alignments, strains)
				alignmentTypes := []string{"core", "disp", "pan"}
				alignmentArray := [][]seqrecord.SeqRecords{cores, disps, alignments}
				for i := 0; i < len(alignmentTypes); i++ {
					name := alignmentTypes[i]
					alns := alignmentArray[i]
					if len(alns) == 0 {
						WARN.Printf("%s, %s alignments has zero record.\n", prefix, name)
					} else {
						j := job{}
						j.alignments = alns
						j.name = name
						j.position = pos
						j.strains = strains
						jobChan <- j
					}
				}
			}
		}
	}()

	for j := range jobChan {
		cmd.RunOne(j.strains, j.alignments, j.position, j.name)
	}
}

func separeteAlignments(alignments []seqrecord.SeqRecords, strains []strain.Strain) (cores, disps []seqrecord.SeqRecords) {
	for _, aln := range alignments {
		genomeSet := make(map[string]bool)
		for _, rec := range aln {
			genomeSet[rec.Genome] = true
		}

		if len(genomeSet) == len(strains) {
			cores = append(cores, aln)
		} else if len(genomeSet) < len(strains) {
			disps = append(disps, aln)
		}
	}

	return
}

func (cmd *cmdCovGenomes) RunOne(strains []strain.Strain, alignments []seqrecord.SeqRecords, pos int, name string) {
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

	ncpu := 1
	done := make(chan bool)
	for i := 0; i < ncpu; i++ {
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
					cc := cov.GenomesCalc(alignments, g, cmd.maxl, pos, covGenomesFunc)
					res := createCovResult(cc, cmd.maxl, pos)
					// Write result to files.
					filePrefix := fmt.Sprintf("%s_%s_%s_pos%d", g.RefAcc(),
						funcType, name, pos)
					filePath := filepath.Join(*cmd.workspace, cmd.covOutBase, s.Path,
						filePrefix+".json")
					if !math.IsNaN(res.VarKs) {
						save2Json(res, filePath)
					} else {
						WARN.Printf("%s: VarKs: NaN\n", filePath)
					}

					if cmd.numBoot > 0 {
						ccChan := cmd.boot(cc, cmd.numBoot)
						resChan := cmd.collectBoot(ccChan, pos, cmd.maxl)
						// Write result to files.
						filePath := filepath.Join(*cmd.workspace, cmd.covOutBase, s.Path,
							filePrefix+"_boot.json.zip")
						f, err := os.Create(filePath)
						if err != nil {
							log.Panicln(err)
						}
						defer f.Close()

						w, _ := flate.NewWriter(f, flate.BestCompression)
						defer w.Close()

						encoder := json.NewEncoder(w)
						for res := range resChan {
							if !math.IsNaN(res.VarKs) {
								if err := encoder.Encode(res); err != nil {
									log.Panicln(err)
								}
							}
						}
					}
				}

			}
			done <- true
		}()
	}

	// Waiting for the job done.
	for i := 0; i < ncpu; i++ {
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

func (cmd *cmdCovGenomes) boot(cc []*cov.Calculators, numBoot int) (ccChan chan []*cov.Calculators) {
	// bootstrapping
	bootJobs := make(chan []int)
	go func() {
		defer close(bootJobs)
		for i := 0; i < numBoot; i++ {
			sample := make([]int, len(cc))
			for j := 0; j < len(sample); j++ {
				sample[j] = rand.Intn(len(cc))
			}
			bootJobs <- sample
		}
	}()

	ccChan = make(chan []*cov.Calculators)
	ncpu := *cmd.ncpu
	done := make(chan bool)
	for i := 0; i < ncpu; i++ {
		go func() {
			for sample := range bootJobs {
				randomCC := []*cov.Calculators{}
				for j := 0; j < len(sample); j++ {
					index := sample[j]
					c := cc[index]
					randomCC = append(randomCC, c)
				}
				ccChan <- randomCC
			}
			done <- true
		}()
	}

	go func() {
		defer close(ccChan)
		for i := 0; i < ncpu; i++ {
			<-done
		}
	}()

	return
}

func (cmd *cmdCovGenomes) collectBoot(ccChan chan []*cov.Calculators, pos, maxl int) (resChan chan CovResult) {
	resChan = make(chan CovResult)
	numWorker := *cmd.ncpu
	done := make(chan bool)
	for i := 0; i < numWorker; i++ {
		go func() {
			for cc := range ccChan {
				res := createCovResult(cc, maxl, pos)
				resChan <- res
			}
			done <- true
		}()
	}

	go func() {
		defer close(resChan)
		for i := 0; i < numWorker; i++ {
			<-done
		}
	}()

	return
}

func createCovResult(cc []*cov.Calculators, maxl, pos int) (res CovResult) {
	biasCorrection := false
	c := cov.NewCalculators(maxl, biasCorrection)
	for i := 0; i < len(cc); i++ {
		c.Append(cc[i])
	}

	// Process and return a cov result.
	res.Ks = c.Ks.Mean.GetResult()
	res.VarKs = c.Ks.Var.GetResult()
	res.N = c.Ks.Mean.GetN()

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
		v := c.TCov.GetResult(index)
		xy := c.TCov.GetMeanXY(index)
		n := c.TCov.GetN(index)
		if !math.IsNaN(v) {
			res.CtIndices = append(res.CtIndices, i)
			res.Ct = append(res.Ct, v)
			res.MeanXY = append(res.MeanXY, xy)
			res.CtN = append(res.CtN, n)
		}

		v = c.SCov.MeanVars[index].Mean.GetResult()
		n = c.SCov.MeanVars[index].Mean.GetN()
		if !math.IsNaN(v) {
			res.CsIndices = append(res.CsIndices, i)
			res.Cs = append(res.Cs, v)
			res.CsN = append(res.CsN, n)
		}

		v = c.MCov.MeanVars[index].Mean.GetResult()
		n = c.MCov.MeanVars[index].Mean.GetN()
		if !math.IsNaN(v) {
			res.CmIndices = append(res.CmIndices, i)
			res.Cm = append(res.Cm, v)
			res.CmN = append(res.CmN, n)
		}

		v = c.RCov.MeanVars[index].Mean.GetResult()
		n = c.RCov.MeanVars[index].Mean.GetN()
		if !math.IsNaN(v) {
			res.CrIndices = append(res.CrIndices, i)
			res.Cr = append(res.Cr, v)
			res.CrN = append(res.CrN, n)
		}
	}
	return
}
