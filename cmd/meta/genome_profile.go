package main

import (
	"github.com/mingzhi/meta/strain"
	"github.com/mingzhi/ncbiftp/taxonomy"
)

// Command to do genome position profiling.
type cmdGenomeProfile struct {
	cmdConfig
}

func (cmd *cmdGenomeProfile) Init() {
	// Parse config file and settings.
	cmd.ParseConfig()
	// Load species map.
	cmd.LoadSpeciesMap()
}

func (cmd *cmdGenomeProfile) Run(args []string) {
	cmd.Init()

	jobs := make(chan strain.Strain)
	go func() {
		defer close(jobs)
		for _, strains := range cmd.speciesMap {
			for _, s := range strains {
				jobs <- s
			}
		}
	}()

	base := cmd.refBase
	gcMap := taxonomy.GeneticCodes()

	done := make(chan bool)
	for i := 0; i < *cmd.ncpu; i++ {
		go func() {
			for s := range jobs {
				s.ProfileGenomes(base, gcMap)
			}
			done <- true
		}()
	}

	for i := 0; i < *cmd.ncpu; i++ {
		<-done
	}
}
