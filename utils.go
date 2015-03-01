package meta

import (
	"encoding/json"
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
