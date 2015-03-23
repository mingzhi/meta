package meta

import (
	"encoding/csv"
	"github.com/mingzhi/ncbiutils"
	"io/ioutil"
	"log"
	"os"
	"path/filepath"
	"strconv"
	"strings"
)

var (
	Info *log.Logger
	Warn *log.Logger
)

// GenerateStrainInfors:
// generates strain informations.
// workspace: output dir
// ref: reference genome folder.
// tax: taxonomy ftp folder.
func GenerateStrainInfors(workspace, ref, tax string) (strains []Strain) {

	records := readSummary(ref)

	strainMap := make(map[string]Strain)
	taxMap := ncbiutils.ReadTaxas(tax)
	uidMap := findPaths(ref)

	for _, r := range records {
		s, found := strainMap[r.ProjectId]

		good := true

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
				good = false
			}

			if _, found := uidMap[s.ProjectId]; found {
				s.Path = uidMap[s.ProjectId]
			} else {
				Warn.Printf("Cound not find path for %s\n", s.Name)
				good = false
			}

		}

		if good {
			g := Genome{}
			g.Accession = r.Accession
			g.Replicon = r.Replicon
			g.Length = r.Length
			s.Genomes = append(s.Genomes, g)

			strainMap[r.ProjectId] = s
		}
	}

	for _, s := range strainMap {
		s.Species = strings.Replace(s.Species, " ", "_", -1)
		strains = append(strains, s)
	}

	return
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
// ref: reference genome folder.
func readSummary(ref string) (records []Record) {
	fileName := "summary.txt"
	filePath := filepath.Join(ref, fileName)
	f, err := os.Open(filePath)
	if err != nil {
		Info.Println(filePath)
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

// findPaths returns a map of strain and path map.
func findPaths(ref string) map[string]string {
	m := make(map[string]string)
	fileInfos, err := ioutil.ReadDir(ref)
	if err != nil {
		log.Panic(err)
	}
	for _, fi := range fileInfos {
		if fi.IsDir() {
			terms := strings.Split(fi.Name(), "_")
			if strings.Contains(terms[len(terms)-1], "uid") {
				uid := strings.Replace(terms[len(terms)-1], "uid", "", -1)
				m[uid] = fi.Name()
			}
		}
	}

	return m
}
