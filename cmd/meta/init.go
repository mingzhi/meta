package main

import (
	"encoding/json"
	"github.com/mingzhi/meta"
	"os"
	"path/filepath"
)

// Create a species_map.json.
// If a prefix is given, generate a json file
// for a list of strains.
type cmdInit struct {
	cmdConfig // embed cmdConfig
}

func (cmd *cmdInit) Run(args []string) {
	// Parse config and settings.
	cmd.ParseConfig()
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

	speciesMap := make(map[string][]meta.Strain)
	for _, s := range strains {
		if s.Species != "" {
			speciesMap[s.Species] = append(speciesMap[s.Species], s)
		}
	}

	ss, found := speciesMap[cmd.prefix]
	if found {
		fileName := filepath.Join(*cmd.workspace, cmd.prefix+"_strains.json")
		w, err := os.Create(fileName)
		if err != nil {
			ERROR.Fatalln(err)
		}
		defer w.Close()

		encoder := json.NewEncoder(w)
		err = encoder.Encode(ss)
		if err != nil {
			ERROR.Fatalln(err)
		}
	} else {
		WARN.Printf("Cannot find strains for %s\n", cmd.prefix)
	}
}

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
