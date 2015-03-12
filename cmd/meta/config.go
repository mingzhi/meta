package main

import (
	"flag"
	"github.com/spf13/viper"
	"runtime"
	"strconv"
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
	cmd.config = fs.String("c", "config", "configure file name.")
	cmd.ncpu = fs.Int("ncpu", runtime.NumCPU(), "number of CPUs for using.")
	return fs
}

// Parse configs.
func (cmd *cmdConfig) ParseConfig() {
	// Register viper for configurations.
	viper.SetConfigName(*cmd.config)
	viper.AddConfigPath(*cmd.workspace)
	viper.ReadInConfig()

	// Read settings.
	cmd.refBase = viper.GetString("Reference_Genome_Directory")
	cmd.taxBase = viper.GetString("Taxonomy_Diretory")
	cmd.strainFileName = viper.GetString("Strain_File_Name")
	cmd.pairedEndReadFile1 = viper.GetString("Paired_End_Reads_1")
	cmd.pairedEndReadFile2 = viper.GetString("Paired_End_Reads_2")
	cmd.bowtieThreadsNum = viper.GetInt("Bowtie_Threads_Num")
	cmd.samOutBase = viper.GetString("Sam_Output_Diretory")

	cmd.orthoOutBase = viper.GetString("Ortholog_Output_Diretory")
	cmd.prefix = viper.GetString("Prefix")

	cmd.covReadsFuncName = viper.GetString("Cov_Reads_Func_Name")
	cmd.covOutBase = viper.GetString("Cov_Output_Directory")
	cmd.maxl = viper.GetInt("Cov_Max_Length")
	positions := viper.GetStringSlice("Positions")
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
}
