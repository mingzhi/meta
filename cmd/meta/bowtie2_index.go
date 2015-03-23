package main

import (
	"bytes"
	"github.com/mingzhi/meta"
	"os"
	"os/exec"
	"path/filepath"
	"time"
)

// Indexing reference genome for efficiently read mapping.

type cmdIndex struct {
	cmdConfig // embed cmdConfig.
}

func (cmd *cmdIndex) Run(args []string) {
	// Parse configure and settings.
	cmd.ParseConfig()
	cmd.LoadSpeciesMap()

	// For each strain,
	// 1. check if genome fasta are already index by bowtie-build,
	// 2. if not, send to job chan for building index.
	jobs := make(chan meta.Strain)
	go func() {
		for _, strains := range cmd.speciesMap {
			for _, strain := range strains {
				if isBowtieIndexExist(strain, cmd.refBase) {
					INFO.Printf("%s has already been indexed!\n", strain.Path)
				} else {
					jobs <- strain
				}
			}
		}

		close(jobs)
	}()

	// Create cmd.ncpu worker for building index,
	// send done signal when the job is done.
	done := make(chan bool)
	for i := 0; i < *cmd.ncpu; i++ {
		go func() {
			for strain := range jobs {
				bowtieBuildIndex(strain, cmd.refBase)
			}
			done <- true
		}()
	}

	for i := 0; i < *cmd.ncpu; i++ {
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
	for _, g := range strain.Genomes {
		acc := meta.FindRefAcc(g.Accession)
		genomeFastaPath := filepath.Join(refDir, strain.Path, acc+".fna")
		genomeIndexBase := filepath.Join(refDir, strain.Path, acc)
		cmd := exec.Command("bowtie2-build", "-f", genomeFastaPath, genomeIndexBase)
		// capture stderr.
		stderr := new(bytes.Buffer)
		cmd.Stderr = stderr
		if err := cmd.Run(); err != nil {
			ERROR.Fatalln(string(stderr.Bytes()))
		}
	}
}
