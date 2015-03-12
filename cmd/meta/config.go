package main

import (
	"flag"
	"github.com/jacobstr/confer"
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
	refBase        string // reference genome folder.
	taxBase        string // taxonomy database folder.
	samOutBase     string // sam output folder.
	strainFileName string // strain file name.

	// For align_reads.
	pairedEndReadFile1 string // paired-end read file 1.
	pairedEndReadFile2 string // paired-end read file 2.
	bowtieThreadsNum   int    // bowtie threads number.

	// For cov calculations.
	positions        []int  // positions in genomic profile to be calculated.
	maxl             int    // max length of correlations.
	covReadsFuncName string // cov calculation function name.
	covOutBase       string // cov output folder.

	// For orthoMCL.
	orthoOutBase string // orthologs and alignment output folder.
	prefix       string // prefix
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
	cmd.refBase = config.GetString("genome.reference")
	cmd.taxBase = config.GetString("genome.taxonomy")
	cmd.prefix = config.GetString("species.name")
	cmd.strainFileName = config.GetString("species.file")
	cmd.pairedEndReadFile1 = config.GetString("reads.paired1")
	cmd.pairedEndReadFile2 = config.GetString("reads.paired2")
	cmd.covOutBase = config.GetString("out.cov")
	cmd.samOutBase = config.GetString("out.sam")
	cmd.orthoOutBase = config.GetString("out.ortho")
	cmd.maxl = config.GetInt("cov.maxl")
	cmd.covReadsFuncName = config.GetString("cov.func")
	cmd.bowtieThreadsNum = config.GetInt("bowtie2.threads")

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
