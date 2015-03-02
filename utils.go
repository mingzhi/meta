package meta

import (
	"encoding/json"
	"github.com/mingzhi/ncbiutils"
	"log"
	"os"
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
