package main

import (
	"encoding/json"
	"flag"
	"github.com/mingzhi/meta"
	"github.com/spf13/viper"
	"log"
	"os"
	"path/filepath"
)

// Initialze project.
// Create a species_map.json.

type cmdInit struct {
	workspace *string // workspace.
	config    *string // configure file name.
	prefix    *string // prefix of the species.

	refBase string // reference genome diretory.
	taxBase string // taxonomy database diretory.
}

func (cmd *cmdInit) Flags(fs *flag.FlagSet) *flag.FlagSet {
	cmd.config = fs.String("c", "config", "configure file name")
	cmd.workspace = fs.String("w", "", "workspace")
	cmd.prefix = fs.String("p", "", "prefix")

	return fs
}

func (cmd *cmdInit) init() {
	viper.SetConfigName(*cmd.config)
	viper.AddConfigPath(*cmd.workspace)
	viper.ReadInConfig()

	cmd.refBase = viper.GetString("Reference_Genome_Directory")
	cmd.taxBase = viper.GetString("Taxonomy_Diretory")

	// register logger for meta package.
	registerLogger()
}

func (cmd *cmdInit) Run(args []string) {
	cmd.init()
	var strains []meta.Strain
	if cmd.isSpeciesMapExist() {
		f, err := os.Open(filepath.Join(*cmd.workspace, "reference_strains.json"))
		if err != nil {
			log.Panic(err)
		}
		defer f.Close()
		decoder := json.NewDecoder(f)
		err = decoder.Decode(&strains)
		if err != nil {
			log.Panic(err)
		}
	} else {
		strains = meta.GenerateStrainInfors(*cmd.workspace, cmd.refBase, cmd.taxBase)
		w, err := os.Create(filepath.Join(*cmd.workspace, "reference_strains.json"))
		if err != nil {
			log.Panic(err)
		}
		defer w.Close()
		encoder := json.NewEncoder(w)
		err = encoder.Encode(strains)
		if err != nil {
			log.Panic(err)
		}
	}

	speciesMap := make(map[string][]meta.Strain)
	for _, s := range strains {
		if s.Species != "" {
			speciesMap[s.Species] = append(speciesMap[s.Species], s)
		}
	}

	ss, found := speciesMap[*cmd.prefix]
	if found {
		fileName := filepath.Join(*cmd.workspace, *cmd.prefix+"_strains.json")
		w, err := os.Create(fileName)
		if err != nil {
			log.Panic(err)
		}
		defer w.Close()

		encoder := json.NewEncoder(w)
		err = encoder.Encode(ss)
		if err != nil {
			log.Panic(err)
		}
	} else {
		log.Printf("Cannot find strains for %s\n", *cmd.prefix)
	}
}

func (cmd *cmdInit) isSpeciesMapExist() (isExist bool) {
	filePath := filepath.Join(*cmd.workspace, "reference_strains.json")
	if fi1, err := os.Stat(filePath); err != nil {
		if os.IsNotExist(err) {
			isExist = false
		} else {
			log.Panic(err)
		}
	} else {
		summaryFile := filepath.Join(cmd.refBase, "summary.txt")
		fi2, err := os.Stat(summaryFile)
		if err != nil {
			if os.IsNotExist(err) {
				log.Fatalln("Cannot find summary.txt in reference genome directory!")
			} else {
				log.Panic(err)
			}
		} else {
			if fi1.ModTime().After(fi2.ModTime()) {
				isExist = true
			} else {
				isExist = false
			}
		}
	}

	return
}
