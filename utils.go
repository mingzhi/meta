package meta

import (
	"encoding/json"
	"github.com/mingzhi/ncbiutils"
	"log"
	"os"
	"strings"
)

func ReadSpeciesMap(fileName string) map[string][]Strain {
	m := make(map[string][]Strain)
	f, err := os.Open(fileName)
	if err != nil {
		log.Panic(err)
	}
	defer f.Close()

	dc := json.NewDecoder(f)
	err = dc.Decode(&m)
	if err != nil {
		log.Panic(err)
	}

	return m
}

func ReadAlignments(fileName string) []ncbiutils.SeqRecords {
	records := []ncbiutils.SeqRecords{}
	f, err := os.Open(fileName)
	if err != nil {
		log.Panic(err)
	}
	defer f.Close()

	dc := json.NewDecoder(f)
	err = dc.Decode(&records)
	if err != nil {
		log.Panic(err)
	}

	return records
}

// Read strain information from a json file.
func ReadStrains(fileName string) (strains []Strain) {
	f, err := os.Open(fileName)
	if err != nil {
		log.Panic(err)
	}
	defer f.Close()

	decoder := json.NewDecoder(f)
	err = decoder.Decode(&strains)
	if err != nil {
		log.Panic(err)
	}

	return
}

// Find and clean reference genome accession.
func FindRefAcc(name string) string {
	return (strings.Split(name, ".")[0])
}
