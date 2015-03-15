package main

import (
	"encoding/json"
	"fmt"
	"github.com/mingzhi/meta"
	"os"
	"path/filepath"
	"sort"
)

// Create a species_map.json.
// If a prefix is given, generate a json file
// for a list of strains.
type cmdInit struct {
	cmdConfig // embed cmdConfig
}

// Run command.
// It requires
// 1. workspace;
// 2. reference genome database folder;
// 3. taxonomy database folder.
// 4. species name.
func (cmd *cmdInit) Run(args []string) {
	// Parse config and settings.
	cmd.ParseConfig()

	// If there is no file containing strain informations,
	// or the one is older than summary.txt,
	// create a new one.
	var strains []meta.Strain
	if cmd.isSpeciesMapExist() {
		f, err := os.Open(filepath.Join(*cmd.workspace, "reference_strains.json"))
		if err != nil {
			ERROR.Fatalln(err)
		}
		defer f.Close()
		decoder := json.NewDecoder(f)
		err = decoder.Decode(&strains)
		if err != nil {
			ERROR.Fatalln(err)
		}
	} else {
		strains = meta.GenerateStrainInfors(*cmd.workspace, cmd.refBase, cmd.taxBase)
		w, err := os.Create(filepath.Join(*cmd.workspace, "reference_strains.json"))
		if err != nil {
			ERROR.Fatalln(err)
		}
		defer w.Close()
		encoder := json.NewEncoder(w)
		err = encoder.Encode(strains)
		if err != nil {
			ERROR.Fatalln(err)
		}
	}

	// Create a specie: []strain map.
	speciesMap := make(map[string][]meta.Strain)
	for _, s := range strains {
		if s.Species != "" {
			speciesMap[s.Species] = append(speciesMap[s.Species], s)
		}
	}

	speciesNames := []string{}
	for name, _ := range speciesMap {
		speciesNames = append(speciesNames, name)
	}

	sort.Strings(speciesNames)
	yamlFilePath := filepath.Join(*cmd.workspace, "reference_species.yaml")
	w, err := os.Create(yamlFilePath)
	if err != nil {
		ERROR.Panicf("Cannot create file %s: %v", yamlFilePath, err)
	}
	defer w.Close()
	for _, name := range speciesNames {
		w.WriteString(fmt.Sprintf("%s:\n", name))
		strains := speciesMap[name]
		for _, s := range strains {
			w.WriteString(fmt.Sprintf(" - %s\n", s.Path))
		}
	}
}

// Check if there exists a file containing strain information.
// If exists, also check the modified time, which should be
// after that of summary.txt.
func (cmd *cmdInit) isSpeciesMapExist() (isExist bool) {
	filePath := filepath.Join(*cmd.workspace, "reference_strains.json")
	if fi1, err := os.Stat(filePath); err != nil {
		if os.IsNotExist(err) {
			isExist = false
		} else {
			ERROR.Fatalln(err)
		}
	} else {
		summaryFile := filepath.Join(cmd.refBase, "summary.txt")
		fi2, err := os.Stat(summaryFile)
		if err != nil {
			if os.IsNotExist(err) {
				ERROR.Fatalln("Cannot find summary.txt in reference genome directory!")
			} else {
				ERROR.Fatalln(err)
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
