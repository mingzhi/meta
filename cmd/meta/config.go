package main

import (
	"encoding/json"
	"flag"
	"github.com/jacobstr/confer"
	"github.com/mingzhi/meta"
	"gopkg.in/yaml.v2"
	"io/ioutil"
	"os"
	"path/filepath"
	"runtime"
	"strconv"
	"strings"
)

// Config to read flags and configure file.
type cmdConfig struct {
	// Flags.
	workspace *string // workspace.
	config    *string // configure file name.
	ncpu      *int    // number of CPUs for using.

	// Data diretory and path.
	refBase    string // reference genome folder.
	taxBase    string // taxonomy database folder.
	samOutBase string // sam output folder.

	// For align_reads.
	pairedEndReadFile1   string // paired-end read file 1.
	pairedEndReadFile2   string // paired-end read file 2.
	bowtieThreadsNum     int    // bowtie threads number.
	maximumMismatchCount int    // maximum mismatch count.

	// For cov calculations.
	positions         []int    // positions in genomic profile to be calculated.
	maxl              int      // max length of correlations.
	covReadsFunctions []string // cov calculation function name.
	covOutBase        string   // cov output folder.

	// For orthoMCL.
	orthoOutBase string // orthologs and alignment output folder.

	// Species strain information.
	// prefix: []strain.
	speciesFile string // species YAML file.
	speciesMap  map[string][]meta.Strain
}

func (cmd *cmdConfig) Flags(fs *flag.FlagSet) *flag.FlagSet {
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
		ERROR.Fatalln(err)
	}
	// Automatic binding.
	config.AutomaticEnv()
	// Data
	cmd.refBase = config.GetString("genome.reference")
	cmd.taxBase = config.GetString("genome.taxonomy")
	// Reads
	cmd.pairedEndReadFile1 = config.GetString("reads.paired1")
	cmd.pairedEndReadFile2 = config.GetString("reads.paired2")
	// Outputs
	cmd.covOutBase = config.GetString("out.cov")
	cmd.samOutBase = config.GetString("out.sam")
	cmd.orthoOutBase = config.GetString("out.ortho")
	// Correlation settings
	cmd.maxl = config.GetInt("cov.maxl")
	cmd.covReadsFunctions = config.GetStringSlice("cov.functions")
	// Bowtie2 settings
	cmd.bowtieThreadsNum = config.GetInt("bowtie2.threads")
	cmd.maximumMismatchCount = config.GetInt("bowtie2.maximum_mismatch_count")
	// Bacterial species informations
	cmd.speciesFile = config.GetString("species.file")
	cmd.LoadSpeciesMap()

	positions := config.GetStringSlice("cov.positions")
	for _, p := range positions {
		pos, err := strconv.Atoi(p)
		if err != nil {
			ERROR.Fatalf("Can not convert %s to integer!", p)
		}
		cmd.positions = append(cmd.positions, pos)
	}

	if len(cmd.positions) == 0 {
		WARN.Println("Use default position: 4!")
	}

	if cmd.bowtieThreadsNum <= 0 {
		cmd.bowtieThreadsNum = 1
	}

	runtime.GOMAXPROCS(*cmd.ncpu)
	registerLogger()
}

// Check if there exists a file containing strain information.
// If exists, also check the modified time, which should be
// after that of summary.txt.
func (cmd *cmdConfig) IsReferenceStrainsExist() (isExist bool) {
	filePath := filepath.Join(*cmd.workspace, "reference_strains.json")
	if fi1, err := os.Stat(filePath); err != nil {
		if os.IsNotExist(err) {
			isExist = false
		} else {
			ERROR.Fatalln(err)
		}
	} else {
		summaryFile := filepath.Join(cmd.refBase, "summary.txt")
		fi2, err := os.Stat(summaryFile)
		if err != nil {
			if os.IsNotExist(err) {
				ERROR.Fatalln("Cannot find summary.txt in reference genome directory!")
			} else {
				ERROR.Fatalln(err)
			}
		} else {
			if fi1.ModTime().After(fi2.ModTime()) {
				isExist = true
			} else {
				isExist = false
			}
		}
	}

	return
}

// Read reference_strains.json.
func (cmd *cmdConfig) ReadReferenceStrains() (strains []meta.Strain) {
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

// Create reference_strains.json
func (cmd *cmdConfig) CreateReferenceStrains() (strains []meta.Strain) {
	strains = meta.GenerateStrainInfors(*cmd.workspace, cmd.refBase, cmd.taxBase)
	w, err := os.Create(filepath.Join(*cmd.workspace, "reference_strains.json"))
	if err != nil {
		ERROR.Fatalln(err)
	}
	defer w.Close()
	encoder := json.NewEncoder(w)
	err = encoder.Encode(strains)
	if err != nil {
		ERROR.Fatalln(err)
	}
	return
}

// Read species file,
// return a map[string][]string,
// in which prefix: []strain.path
func (cmd *cmdConfig) ReadSpeciesFile() map[string][]string {
	filePath := filepath.Join(*cmd.workspace, cmd.speciesFile)
	f, err := os.Open(filePath)
	if err != nil {
		ERROR.Fatalf("Cannot open %s: %v\n", filePath, err)
	}
	defer f.Close()

	data, err := ioutil.ReadAll(f)
	if err != nil {
		ERROR.Fatalf("Cannot read file %s: %v\n", filePath, err)
	}

	m := make(map[string][]string)
	err = yaml.Unmarshal(data, &m)
	if err != nil {
		ERROR.Fatalf("Cannot unmarshal %s: %v\n", filePath, err)
	}

	return m
}

// Load species map.
func (cmd *cmdConfig) LoadSpeciesMap() {
	// If reference_strains exists, read it,
	// else generate it.
	var strains []meta.Strain
	if cmd.IsReferenceStrainsExist() {
		strains = cmd.ReadReferenceStrains()
	} else {
		strains = cmd.CreateReferenceStrains()
	}

	// Make a strain map,
	// strain.path: strain
	strainMap := make(map[string]meta.Strain)
	for _, strain := range strains {
		strainMap[strain.Path] = strain
	}

	// Read species input yaml file.
	inputMap := cmd.ReadSpeciesFile()

	// Load species map.
	cmd.speciesMap = make(map[string][]meta.Strain)
	for prefix, strainPaths := range inputMap {
		for _, strainPath := range strainPaths {
			strain, found := strainMap[strainPath]
			if found {
				cmd.speciesMap[prefix] = append(cmd.speciesMap[prefix], strain)
			}
		}
	}
}
