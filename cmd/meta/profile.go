package main

import (
	"flag"
	"github.com/mingzhi/meta"
	"log"
	"path/filepath"
)

type cmdProfile struct {
	workspace *string
	prefix    *string
	ref       *string
}

func (cmd *cmdProfile) Flags(fs *flag.FlagSet) *flag.FlagSet {
	cmd.workspace = fs.String("w", "", "workspace")
	cmd.ref = fs.String("r", "", "reference genome diretory")
	cmd.prefix = fs.String("p", "", "prefix")

	return fs
}

func (cmd *cmdProfile) Run(args []string) {
	registerLogger()
	speciesMapFileName := filepath.Join(*cmd.workspace, "species_map.json")
	speciesMap := meta.ReadSpeciesMap(speciesMapFileName)
	strains, found := speciesMap[*cmd.prefix]
	if !found {
		log.Fatalf("Can not find strain information for %s from species_map.json\n", *cmd.prefix)
	}

	meta.GenomePosProfiling(strains, *cmd.ref)
}
