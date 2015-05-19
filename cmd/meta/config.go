package main

import (
	"encoding/json"
	"flag"
	"github.com/jacobstr/confer"
	"github.com/mingzhi/meta/strain"
	"gopkg.in/yaml.v2"
	"io/ioutil"
	"os"
	"path/filepath"
	"runtime"
	"strconv"
	"strings"
)

type fitControl struct {
	name       string
	start, end int
}

// Config to read flags and configure file.
type cmdConfig struct {
	// Flags.
	workspace *string // workspace.
	config    *string // configure file name.
	ncpu      *int    // number of CPUs for using.

	// Data diretory and path.
	refBase string // reference genome folder.
	taxBase string // taxonomy database folder.
	repBase string // genome report database.

	// Output folders.
	samOutBase   string // sam output folder.
	covOutBase   string // cov output folder.
	orthoOutBase string // orthologs and alignment output folder.
	plotOutBase  string // plot output folder.
	fitOutBase   string // fitting output folder.

	// For align_reads.
	// Now we only support paired-end reads.
	pairedEndReadFile1 string // paired-end read file 1.
	pairedEndReadFile2 string // paired-end read file 2.

	// Bowtie2 options.
	bowtieOptions []string // bowtie2 options.

	// For cov calculations.
	positions     []int    // positions in genomic profile to be calculated.
	maxl          int      // max length of correlations.
	covReadsFuncs []string // cov calculation function name.

	// Species strain information.
	speciesFile string                     // species YAML file.
	speciesMap  map[string][]strain.Strain // species: []strain map.

	// Fit parameters.
	fitControls []fitControl

	// bootstrapping parameters.
	numBoot int // number of bootstrapping
}

// Implement command package interface.
func (cmd *cmdConfig) Flags(fs *flag.FlagSet) *flag.FlagSet {
	// define subcommand's flags.
	cmd.workspace = fs.String("w", "", "workspace.")
	cmd.config = fs.String("c", "config.yaml", "configure files in YAML format, which are separeted by comma.")
	cmd.ncpu = fs.Int("ncpu", runtime.NumCPU(), "number of CPUs for using.")
	return fs
}

// Parse configs.
func (cmd *cmdConfig) ParseConfig() {
	// Use confer package to parse configure files.
	config := confer.NewConfig()
	// Set root path, which contains configure files.
	config.SetRootPath(*cmd.workspace)
	// Read configure files.
	configPaths := strings.Split(*cmd.config, ",")
	if err := config.ReadPaths(configPaths...); err != nil {
		ERROR.Panicln(err)
	}
	// Automatic binding.
	config.AutomaticEnv()
	// Parse data directory and path.
	cmd.refBase = config.GetString("genome.reference")
	cmd.taxBase = config.GetString("genome.taxonomy")
	cmd.repBase = config.GetString("genome.reports")
	// Parse file names of paried-end reads.
	cmd.pairedEndReadFile1 = config.GetString("reads.paired1")
	cmd.pairedEndReadFile2 = config.GetString("reads.paired2")
	// Parse output folders.
	cmd.covOutBase = config.GetString("out.cov")
	cmd.samOutBase = config.GetString("out.sam")
	cmd.orthoOutBase = config.GetString("out.ortho")
	cmd.plotOutBase = config.GetString("out.plot")
	// Parse options for covariance calculation.
	cmd.maxl = config.GetInt("cov.maxl")
	cmd.covReadsFuncs = config.GetStringSlice("cov.functions")
	// Parse positions to be calculated.
	positions := config.GetStringSlice("cov.positions")
	for _, p := range positions {
		pos, err := strconv.Atoi(p)
		if err != nil {
			ERROR.Panicf("Can not convert %s to integer!", p)
		}
		cmd.positions = append(cmd.positions, pos)
	}
	// Parse bowtie2 options.
	cmd.bowtieOptions = config.GetStringSlice("bowtie2.options")
	// Parse the name of file storing bacterial strain information.
	cmd.speciesFile = config.GetString("species.file")

	// Fit range.
	expFitControl := fitControl{}
	expFitControl.start = config.GetInt("fit.exp.start")
	expFitControl.end = config.GetInt("fit.exp.end")
	expFitControl.name = "exp"
	cmd.fitControls = append(cmd.fitControls, expFitControl)
	hyperFitControl := fitControl{}
	hyperFitControl.start = config.GetInt("fit.hyper.start")
	hyperFitControl.end = config.GetInt("fit.hyper.end")
	hyperFitControl.name = "hyper"
	cmd.fitControls = append(cmd.fitControls, hyperFitControl)

	// Bootstrapping
	cmd.numBoot = config.GetInt("bootstrapping.number")

	runtime.GOMAXPROCS(*cmd.ncpu)
}

// Read reference_strains.json.
func (cmd *cmdConfig) ReadReferenceStrains() (strains []strain.Strain) {
	f, err := os.Open(filepath.Join(*cmd.workspace, "reference_strains.json"))
	if err != nil {
		ERROR.Fatalln(err)
	}
	defer f.Close()
	decoder := json.NewDecoder(f)
	err = decoder.Decode(&strains)
	if err != nil {
		ERROR.Fatalln(err)
	}
	return
}

// Read species file in yaml file format.
// return a map[string][]string,
// in which prefix: []strain.path
func (cmd *cmdConfig) ReadSpeciesFile() map[string][]string {
	filePath := filepath.Join(*cmd.workspace, cmd.speciesFile)
	f, err := os.Open(filePath)
	if err != nil {
		ERROR.Panicln("Cannot open %s: %v\n", filePath, err)
	}
	defer f.Close()

	data, err := ioutil.ReadAll(f)
	if err != nil {
		ERROR.Panicln("Cannot read file %s: %v\n", filePath, err)
	}

	m := make(map[string][]string)
	err = yaml.Unmarshal(data, &m)
	if err != nil {
		ERROR.Panicln("Cannot unmarshal %s: %v\n", filePath, err)
	}

	return m
}

// Load species strains map.
func (cmd *cmdConfig) LoadSpeciesMap() {
	// If reference_strains exists, read it.
	var strains []strain.Strain
	if isReferenceStrainsExists(*cmd.workspace, cmd.repBase) {
		strains = cmd.ReadReferenceStrains()
	} else {
		ERROR.Fatalln("Can not find reference_strains.json, please run meta init first!")
	}

	// Make a strain map,
	// strain.path: strain
	strainMap := make(map[string]strain.Strain)
	for _, strain := range strains {
		strainMap[strain.Path] = strain
	}

	// Read species input yaml file.
	inputMap := cmd.ReadSpeciesFile()

	// Load species map.
	cmd.speciesMap = make(map[string][]strain.Strain)
	for prefix, strainPaths := range inputMap {
		for _, strainPath := range strainPaths {
			strain, found := strainMap[strainPath]
			if found {
				if cmd.checkStrain(strain) {
					cmd.speciesMap[prefix] = append(cmd.speciesMap[prefix], strain)
				}
			} else {
				WARN.Printf("Cannot find %s\n", strainPath)
			}
		}
	}
}

// check if the strain has complete genome sequences.
func (cmd *cmdConfig) checkStrain(s strain.Strain) bool {
	// for each genome, check if its .fna file exists.
	for _, g := range s.Genomes {
		fileName := g.RefAcc() + ".fna"
		filePath := filepath.Join(cmd.refBase, s.Path, fileName)
		if _, err := os.Stat(filePath); err != nil {
			WARN.Println(err)
			return false
		}
	}

	return true
}
