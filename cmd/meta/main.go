package main

import (
	"flag"
	"github.com/mingzhi/meta"
	"github.com/rakyll/command"
	"log"
	"os"
	"path/filepath"
	"runtime"
	"strings"
)

var (
	DefaultMaxProcs = runtime.NumCPU()
)

func main() {
	runtime.GOMAXPROCS(DefaultMaxProcs)
	command.On("init", "initialize", &initCmd{}, []string{})
	command.On("makeudb", "make usearch db", &makeUdbCmd{}, []string{})

	command.ParseAndRun()
}

func registerLogger() {
	meta.Info = log.New(os.Stdout, "INFO: ", log.Ldate|log.Ltime|log.Lshortfile)
	meta.Warn = log.New(os.Stdout, "WARN: ", log.Ldate|log.Ltime|log.Lshortfile)
}

type initCmd struct {
	workspace *string
	ref       *string
	tax       *string
}

func (cmd *initCmd) Flags(fs *flag.FlagSet) *flag.FlagSet {
	cmd.workspace = fs.String("w", "", "workspace")
	cmd.ref = fs.String("r", "", "reference genome diretory")
	cmd.tax = fs.String("t", "", "taxonomy ftp dump folder")

	return fs
}

func (cmd *initCmd) Run(args []string) {
	registerLogger()
	meta.GenerateSpeciesMap(*cmd.workspace, *cmd.ref, *cmd.tax)
}

type makeUdbCmd struct {
	workspace *string
	ref       *string
	prefix    *string
}

func (cmd *makeUdbCmd) Flags(fs *flag.FlagSet) *flag.FlagSet {
	cmd.workspace = fs.String("w", "", "workspace")
	cmd.ref = fs.String("r", "", "reference genome diretory")
	cmd.prefix = fs.String("p", "", "prefix")

	return fs
}

func (cmd *makeUdbCmd) Run(args []string) {
	registerLogger()
	speciesMapFileName := filepath.Join(*cmd.workspace, "species_map.json")
	speciesMap := meta.ReadSpeciesMap(speciesMapFileName)
	strains, found := speciesMap[*cmd.prefix]
	if !found {
		log.Fatalf("Can not find strain information for %s from species_map.json\n", *cmd.prefix)
	}

	for _, s := range strains {
		for _, g := range s.Genomes {
			acc := cleanAccession(g.Accession)
			fileName := filepath.Join(*cmd.ref, s.Path, acc+".faa")
			meta.UsearchMakeUDB(fileName)
		}
	}
}

func cleanAccession(name string) string {
	return strings.Split(name, ".")[0]
}
