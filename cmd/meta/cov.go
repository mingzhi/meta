package main

// Caclulate sequence correlations for each reference genome,
// by mapping reads to it.

import (
	"bufio"
	"encoding/json"
	"flag"
	"fmt"
	"github.com/mingzhi/meta"
	"github.com/spf13/viper"
	"io"
	"log"
	"math"
	"math/rand"
	"os"
	"path/filepath"
	"runtime"
	"strings"
)

type CovResult struct {
	Ks, VarKs float64
	Ct        []float64
	CtIndices []int
	CtN       []int
	N         int
}

type covGenomeFunc func(records meta.SamRecords, genome meta.Genome, maxl int, pos int) (*meta.KsCalculator, *meta.CovCalculator)

type cmdCov struct {
	workspace *string // workspace.
	config    *string // configure file name.
	prefix    *string // prefix of bacterial species.
	ncpu      *int    // number of CPUs for using.

	samout      string // sam output diretory.
	outdir      string
	ref         string        // reference genome database.
	maxl        int           // max length of correlation.
	bootstrap   int           // number of bootstraps.
	covFunc     covGenomeFunc // cov calculate function.
	covFuncName string        // cov calculate function name.
	strainFile  string        // strain file name.

	pos int // position in a codon to calculate.
}

func (cmd *cmdCov) Flags(fs *flag.FlagSet) *flag.FlagSet {
	cmd.config = fs.String("c", "config", "configure file name")
	cmd.workspace = fs.String("w", "", "workspace")
	cmd.prefix = fs.String("p", "", "prefix")
	cmd.ncpu = fs.Int("cpu", runtime.NumCPU(), "Number of CPUs for using")

	return fs
}

func (cmd *cmdCov) init() {
	viper.SetConfigName(*cmd.config)
	viper.AddConfigPath(*cmd.workspace)
	viper.ReadInConfig()

	cmd.ref = viper.GetString("reference")
	cmd.samout = viper.GetString("samout")
	cmd.outdir = viper.GetString("outdir")
	cmd.maxl = viper.GetInt("maxl")
	cmd.bootstrap = viper.GetInt("bootstrap")

	cmd.covFuncName = viper.GetString("CovGenomeFunc")
	switch cmd.covFuncName {
	case "CovReads":
		cmd.covFunc = meta.CovReads
	default:
		cmd.covFunc = meta.CovGenome
		cmd.covFuncName = "CovGenome"
	}

	if *cmd.prefix == "" {
		cmd.strainFile = viper.GetString("strains")
	}

	runtime.GOMAXPROCS(*cmd.ncpu)
}

func (cmd *cmdCov) Run(args []string) {
	cmd.init()
	registerLogger()
	// Load species map.
	speciesMapFileName := filepath.Join(*cmd.workspace, "species_map.json")
	speciesMap := meta.ReadSpeciesMap(speciesMapFileName)
	// Find the wanted species,
	// stop if it could not find.
	strains := cmd.getStrains(speciesMap)
	if len(strains) == 0 {
		log.Fatalf("Did not find any strain for %s from species_map.json\n",
			*cmd.prefix)
	}

	// Make position profiles.
	meta.GenomePosProfiling(strains, cmd.ref)

	positions := []int{0, 1, 2, 3, 4} // positions to be calculated.

	for _, s := range strains {
		for _, genome := range s.Genomes {
			// We now only looks at chromosome genome.
			if strings.Contains(genome.Replicon, "chromosome") {
				// Read read records from a bam file.
				bamFileName := s.Path + ".align.bam"
				bamFilePath := filepath.Join(cmd.samout, s.Path, bamFileName)
				// Check if the bam file exists.
				if _, err := os.Stat(bamFilePath); err != nil {
					if os.IsNotExist(err) {
						WARN.Printf("%s does not exist!\n", bamFilePath)
						continue
					} else {
						log.Panicln(err)
					}
				}
				_, records := meta.ReadBamFile(bamFilePath)
				if len(records) == 0 {
					WARN.Printf("%s\t%s zero records\n", s.Path, genome.Accession)
					continue
				}

				// Read postion profile for the genome.
				posFileName := findRefAcc(genome.Accession) + ".pos"
				posFilePath := filepath.Join(cmd.ref, s.Path, posFileName)
				genome.PosProfile = meta.ReadPosProfile(posFilePath)

				// Read sequence for the genome.
				fnaFileName := findRefAcc(genome.Accession) + ".fna"
				fnaFilePath := filepath.Join(cmd.ref, s.Path, fnaFileName)
				genome.Seq = meta.ReadFasta(fnaFilePath).Seq

				for _, pos := range positions {
					cr := cmd.cov(records, genome, pos)

					filePrefix := fmt.Sprintf("%s_%s_pos%d", s.Path, cmd.covFuncName, pos)

					if !math.IsNaN(cr.VarKs) {
						save2Json(cr, filepath.Join(*cmd.workspace, cmd.outdir, filePrefix+".json"))
						save2tsv(cr, filepath.Join(*cmd.workspace, cmd.outdir, filePrefix+".tsv"))
					}

					// Bootstrap
					if cmd.bootstrap > 0 {
						crChan := cmd.boostrapping(records, genome, pos)
						crs := []CovResult{}
						for cr := range crChan {
							if !math.IsNaN(cr.VarKs) {
								crs = append(crs, cr)
							}
						}
						fileName := filePrefix + "_bootstrap.json"
						filePath := filepath.Join(*cmd.workspace, cmd.outdir, fileName)
						f, err := os.Create(filePath)
						if err != nil {
							log.Panic(err)
						}
						defer f.Close()

						e := json.NewEncoder(f)
						err = e.Encode(crs)
						if err != nil {
							log.Panic(err)
						}
					}
				}
			}
		}
	}
}

// Boostrapping.
func (cmd *cmdCov) boostrapping(records meta.SamRecords, genome meta.Genome, pos int) chan CovResult {
	ncpu := runtime.GOMAXPROCS(0)
	c := make(chan CovResult)
	jobs := make(chan int)
	go func() {
		for i := 0; i < cmd.bootstrap; i++ {
			jobs <- i
		}
		close(jobs)
	}()

	done := make(chan int)
	for i := 0; i < ncpu; i++ {
		go func() {
			for j := range jobs {
				rs := sample(records)
				cr := cmd.cov(rs, genome, pos)
				c <- cr
				done <- j
			}
		}()
	}

	// wait and close result channel.
	go func() {
		for i := 0; i < cmd.bootstrap; i++ {
			<-done
		}
		close(c)
	}()
	return c
}

// Calculate covariance for records.
func (cmd *cmdCov) cov(records meta.SamRecords, genome meta.Genome, pos int) CovResult {
	kc, cc := cmd.covFunc(records, genome, cmd.maxl, pos)

	cr := CovResult{}
	cr.Ks = kc.Mean.GetResult()
	cr.VarKs = kc.Var.GetResult()
	cr.N = kc.Mean.GetN()
	var step, length int
	if pos == 0 {
		step = 1
		length = cmd.maxl
	} else {
		step = 3
		length = cmd.maxl / 3
	}
	for i := 1; i < length; i++ {
		v := cc.GetResult(step * i)
		n := cc.GetN(step * i)
		if !math.IsNaN(v) {
			cr.CtIndices = append(cr.CtIndices, i)
			cr.Ct = append(cr.Ct, v)
			cr.CtN = append(cr.CtN, n)
		}
	}

	return cr
}

// Obtain strains.
func (cmd *cmdCov) getStrains(m map[string][]meta.Strain) (strains []meta.Strain) {
	sM := make(map[string]meta.Strain)
	for _, strains := range m {
		for _, s := range strains {
			sM[s.Path] = s
		}
	}

	if cmd.strainFile != "" {
		f, err := os.Open(cmd.strainFile)
		if err != nil {
			log.Fatalf("Cannot open file: %s\n", cmd.strainFile)
		}
		defer f.Close()

		r := bufio.NewReader(f)
		prefixes := []string{}
		for {
			line, err := r.ReadString('\n')
			if err != nil {
				if err != io.EOF {
					log.Panic(err)
				}
				break
			}

			prefix := strings.TrimSpace(line)
			prefixes = append(prefixes, prefix)
		}

		for _, p := range prefixes {
			s, found := sM[p]
			if found {
				strains = append(strains, s)
			}
		}
	} else {
		strains, _ = m[*cmd.prefix]
		s, found := sM[*cmd.prefix]
		if found {
			strains = append(strains, s)
		}
	}

	return
}

// Sample with replacement.
func sample(records meta.SamRecords) meta.SamRecords {
	newRecords := meta.SamRecords{}
	for i := 0; i < len(records); i++ {
		r := rand.Intn(len(records))
		newRecords = append(newRecords, records[r])
	}
	return newRecords
}

func save2Json(cr CovResult, fileName string) {
	f, err := os.Create(fileName)
	if err != nil {
		log.Panic(err)
	}
	defer f.Close()

	ec := json.NewEncoder(f)
	if err := ec.Encode(cr); err != nil {
		log.Panic(err)
	}
}

func save2tsv(cr CovResult, fileName string) {
	f, err := os.Create(fileName)
	if err != nil {
		log.Panic(err)
	}
	defer f.Close()

	f.WriteString(fmt.Sprintf("#Ks=%f;VarKs=%f;N=%d\n", cr.Ks, cr.VarKs, cr.N))
	for i, index := range cr.CtIndices {
		v := cr.Ct[i]
		n := cr.CtN[i]
		f.WriteString(fmt.Sprintf("%d\t%f\t%d\n", index, v, n))
	}
}
