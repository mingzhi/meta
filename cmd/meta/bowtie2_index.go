package main

import (
	"bytes"
	"flag"
	"github.com/mingzhi/meta"
	"github.com/spf13/viper"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"strings"
	"time"
)

// Indexing reference genome for efficiently read mapping.

type cmdIndex struct {
	workspace *string // workspace.
	config    *string // configure file name.

	refDir         string // reference genome folder.
	strainFileName string // strain file name.
	ncpu           int    // number of CPUs for using.
}

// Register flags.
func (cmd *cmdIndex) Flags(fs *flag.FlagSet) *flag.FlagSet {
	cmd.workspace = fs.String("w", "", "workspace")
	cmd.config = fs.String("c", "config", "configure file name")
	return fs
}

// Initialize.
func (cmd *cmdIndex) Init() {
	// Register viper for configurations.
	viper.SetConfigName(*cmd.config)
	viper.AddConfigPath(*cmd.workspace)
	viper.ReadInConfig()

	// Read settings.
	cmd.refDir = viper.GetString("Reference_Genome_Directory")
	cmd.strainFileName = viper.GetString("Strain_File_Name")
	cmd.ncpu = runtime.GOMAXPROCS(0)
}

func (cmd *cmdIndex) Run(args []string) {
	cmd.Init()
	// Read strain information.
	strainFilePath := filepath.Join(*cmd.workspace, cmd.strainFileName)
	strains := meta.ReadStrains(strainFilePath)
	// For each strain,
	// 1. check if genome fasta are already index by bowtie-build,
	// 2. if not, send to job chan for building index.
	jobs := make(chan meta.Strain)
	go func() {
		for _, strain := range strains {
			if isBowtieIndexExist(strain, cmd.refDir) {
				INFO.Printf("%s has already been indexed!\n", strain.Path)
			} else {
				jobs <- strain
			}
		}
		close(jobs)
	}()

	// Create cmd.ncpu worker for building index,
	// send done signal when the job is done.
	done := make(chan bool)
	for i := 0; i < cmd.ncpu; i++ {
		go func() {
			for strain := range jobs {
				bowtieBuildIndex(strain, cmd.refDir)
			}
			done <- true
		}()
	}

	for i := 0; i < cmd.ncpu; i++ {
		<-done
	}
}

func isBowtieIndexExist(strain meta.Strain, refDir string) (isExist bool) {
	genomeIndexFileName := strain.Path + ".1.bt2"
	genomeIndexFilePath := filepath.Join(refDir, strain.Path,
		genomeIndexFileName)
	// Check if bowtie index file exist.
	if fileInfo, err := os.Stat(genomeIndexFilePath); err != nil {
		if os.IsNotExist(err) {
			isExist = false
		} else {
			ERROR.Fatalln(err)
		}
	} else {
		// Get modified time for index,
		// and compare it to the max modifed time
		// of fasta files.
		indexTime := fileInfo.ModTime()
		var maxFastaTime time.Time
		for _, g := range strain.Genomes {
			acc := meta.FindRefAcc(g.Accession)
			genomeFastaName := acc + ".fna"
			genomeFastaPath := filepath.Join(refDir, strain.Path,
				genomeFastaName)
			if fi, err := os.Stat(genomeFastaPath); err != nil {
				ERROR.Fatalln(err)
			} else {
				if maxFastaTime.Before(fi.ModTime()) {
					maxFastaTime = fi.ModTime()
				}
			}
		}

		if indexTime.After(maxFastaTime) {
			isExist = true
		} else {
			isExist = false
		}
	}

	return
}

// Build bowtie2 index.
func bowtieBuildIndex(strain meta.Strain, refDir string) {
	genomeFastaPaths := []string{}
	for _, g := range strain.Genomes {
		acc := meta.FindRefAcc(g.Accession)
		genomeFastaName := acc + ".fna"
		genomeFastaPath := filepath.Join(refDir, strain.Path, genomeFastaName)
		genomeFastaPaths = append(genomeFastaPaths, genomeFastaPath)
	}

	genomeFastasPath := strings.Join(genomeFastaPaths, ",")
	genomeIndexBase := filepath.Join(refDir, strain.Path, strain.Path)

	cmd := exec.Command("bowtie2-build", "-f", genomeFastasPath, genomeIndexBase)
	// capture stderr.
	stderr := new(bytes.Buffer)
	cmd.Stderr = stderr
	if err := cmd.Run(); err != nil {
		ERROR.Fatalln(string(stderr.Bytes()))
	}
}
