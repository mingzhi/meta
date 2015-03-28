package main

import (
	"bytes"
	"github.com/mingzhi/meta/strain"
	"os/exec"
	"path/filepath"
)

// Indexing reference genome for efficiently read mapping.

type cmdIndex struct {
	cmdConfig // embed cmdConfig.
}

func (cmd *cmdIndex) Run(args []string) {
	// Parse configure and settings.
	cmd.ParseConfig()
	cmd.LoadSpeciesMap()

	jobs := make(chan strain.Strain)
	go func() {
		for _, strains := range cmd.speciesMap {
			for _, strain := range strains {
				jobs <- strain
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

// Build bowtie2 index.
func bowtieBuildIndex(strain strain.Strain, refDir string) {
	for _, g := range strain.Genomes {
		genomeFastaPath := filepath.Join(refDir, strain.Path, g.RefAcc()+".fna")
		genomeIndexBase := filepath.Join(refDir, strain.Path, g.RefAcc())
		cmd := exec.Command("bowtie2-build", "-f", genomeFastaPath, genomeIndexBase)
		// capture stderr.
		stderr := new(bytes.Buffer)
		cmd.Stderr = stderr
		if err := cmd.Run(); err != nil {
			ERROR.Println(string(stderr.Bytes()))
			ERROR.Panic(err)
		}
	}
}
