package main

import (
	"flag"
	"github.com/mingzhi/meta"
	"github.com/spf13/viper"
)

// Initialze project.
// Create a species_map.json.

type cmdInit struct {
	workspace *string // workspace.
	config    *string // configure file name.

	ref string // reference genome database.
	tax string // taxonomy database folder.
}

func (cmd *cmdInit) Flags(fs *flag.FlagSet) *flag.FlagSet {
	cmd.config = fs.String("c", "config", "configure file name")
	cmd.workspace = fs.String("w", "", "workspace")

	return fs
}

func (cmd *cmdInit) init() {
	viper.SetConfigName(*cmd.config)
	viper.AddConfigPath(*cmd.workspace)
	viper.ReadInConfig()

	cmd.ref = viper.GetString("reference")
	cmd.tax = viper.GetString("taxonomy")

	// register logger for meta package.
	registerLogger()
}

func (cmd *cmdInit) Run(args []string) {
	cmd.init()
	meta.GenerateSpeciesMap(*cmd.workspace, cmd.ref, cmd.tax)
}
