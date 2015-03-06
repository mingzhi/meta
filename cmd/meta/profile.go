package main

import (
	"flag"
	"github.com/mingzhi/meta"
	"github.com/spf13/viper"
	"log"
	"path/filepath"
	"runtime"
)

type cmdProfile struct {
	workspace *string // workspace.
	config    *string // configure file name.
	prefix    *string // prefix of bacterial species.
	ncpu      *int    // number of CPUs for using.

	ref string // reference genome database.
}

func (cmd *cmdProfile) Flags(fs *flag.FlagSet) *flag.FlagSet {
	cmd.config = fs.String("c", "config", "configure file name")
	cmd.workspace = fs.String("w", "", "workspace")
	cmd.prefix = fs.String("p", "", "prefix")
	cmd.ncpu = fs.Int("cpu", runtime.NumCPU(), "Number of CPUs for using")

	return fs
}

func (cmd *cmdProfile) init() {
	viper.SetConfigName(*cmd.config)
	viper.AddConfigPath(*cmd.workspace)
	viper.ReadInConfig()

	cmd.ref = viper.GetString("reference")
}

func (cmd *cmdProfile) Run(args []string) {
	registerLogger()
	speciesMapFileName := filepath.Join(*cmd.workspace, "species_map.json")
	speciesMap := meta.ReadSpeciesMap(speciesMapFileName)
	strains, found := speciesMap[*cmd.prefix]
	if !found {
		log.Fatalf("Can not find strain information for %s from species_map.json\n", *cmd.prefix)
	}

	meta.GenomePosProfiling(strains, cmd.ref)
}
