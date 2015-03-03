// This programe generate a species map {speciesName: []Strains}.
package main

import (
	"encoding/csv"
	"encoding/json"
	"flag"
	"github.com/mingzhi/ncbiutils"
	"log"
	"os"
	"path/filepath"
	"strconv"
)

type Strain struct {
	Name        string   // taxonomy name.
	TaxId       string   // taxonomy Id.
	ProjectId   string   // project Id.
	Genomes     []Genome // genomes.
	Path        string   // path in referenc genome diretory.
	GeneticCode string   // genetic code id.
	Species     string   // species name
}

type Genome struct {
	Accession string
	Replicon  string
	Length    int
}

var (
	workspace, ref, tax string
)

var (
	Info *log.Logger
	Warn *log.Logger
)

func init() {
	flag.StringVar(&workspace, "w", "", "workspace")
	flag.StringVar(&ref, "r", "", "reference genome diretory")
	flag.StringVar(&tax, "t", "", "taxonomy ftp diretory")
	Info = log.New(os.Stdout, "INFO: ", log.Ldate|log.Ltime|log.Lshortfile)
	Warn = log.New(os.Stdout, "WARN: ", log.Ldate|log.Ltime|log.Lshortfile)
}

func main() {
	flag.Parse()

	records := readSummary()

	strainMap := make(map[string]Strain)
	taxMap := ncbiutils.ReadTaxas(tax)

	for _, r := range records {
		s, found := strainMap[r.ProjectId]

		if !found {
			s = Strain{}
			s.ProjectId = r.ProjectId
			s.Name = r.TaxName
			s.TaxId = r.TaxId
			if _, found := taxMap[s.TaxId]; found {
				s.GeneticCode = findGeneticCode(s.TaxId, taxMap)
				s.Species = findSpecieName(s.TaxId, taxMap)
			} else {
				Warn.Printf("Could not find taxonomy for %s, %s\n", s.TaxId, s.Name)
			}
		}

		g := Genome{}
		g.Accession = r.Accession
		g.Replicon = r.Replicon
		g.Length = r.Length
		s.Genomes = append(s.Genomes, g)

		strainMap[r.ProjectId] = s
	}

	speciesMap := make(map[string][]Strain)
	for _, s := range strainMap {
		if s.Species != "" {
			speciesMap[s.Species] = append(speciesMap[s.Species], s)
		}
	}

	// save to a json file.
	outFileName := "species_map.json"
	outFilePath := filepath.Join(workspace, outFileName)
	w, err := os.Create(outFilePath)
	if err != nil {
		log.Panic(err)
	}
	defer w.Close()

	ec := json.NewEncoder(w)
	err = ec.Encode(speciesMap)
	if err != nil {
		log.Panic(err)
	}

	Info.Printf("Total %d species.\n", len(speciesMap))
}

type Record struct {
	Accession  string
	GenbankAcc string
	Length     int
	TaxId      string
	ProjectId  string
	TaxName    string
	Replicon   string
}

// Read summary.txt in Bacteria ftp folder.
func readSummary() (records []Record) {
	fileName := "summary.txt"
	filePath := filepath.Join(ref, fileName)
	f, err := os.Open(filePath)
	if err != nil {
		log.Fatalln("Could not find summary.txt!")
	}
	defer f.Close()

	csvReader := csv.NewReader(f)
	csvReader.Comma = '\t'
	rows, err := csvReader.ReadAll()
	// First line is the header:
	// Accession	GenbankAcc	Length	Taxid	ProjectID	TaxName	Replicon	Create Date	Update Date
	for i := 1; i < len(rows); i++ {
		row := rows[i]
		r := Record{}
		r.Accession = row[0]
		r.GenbankAcc = row[1]
		r.Length, _ = strconv.Atoi(row[2])
		r.ProjectId = row[4]
		r.TaxId = row[3]
		r.TaxName = row[5]
		r.Replicon = row[6]
		records = append(records, r)
	}

	return
}

// Given the taxonomy Id,
// search the species name.
func findSpecieName(id string, m map[string]ncbiutils.Taxa) string {
	t := m[id]
	var name string
	for {
		if t.Rank == "species" {
			name = t.Name
			break
		}
		t = m[t.Parent]
	}

	return name
}

// Given the taxonomy Id,
// return the genetic code Id.
func findGeneticCode(id string, m map[string]ncbiutils.Taxa) string {
	return m[id].GeneticCode.Id
}
