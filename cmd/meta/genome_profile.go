package main

import (
	"github.com/mingzhi/meta"
)

// Command to do genome position profiling.
type cmdGenomeProfile struct {
	cmdConfig
}

func (cmd *cmdGenomeProfile) Run(args []string) {
	cmd.ParseConfig()
	cmd.LoadSpeciesMap()
	for _, strains := range cmd.speciesMap {
		meta.GenomePosProfiling(strains, cmd.refBase)
	}
}
