package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"github.com/mingzhi/meta/genome"
	"github.com/mingzhi/meta/strain"
	"github.com/mingzhi/ncbiftp/genomes/reports"
	"github.com/mingzhi/ncbiftp/taxonomy"
	"io/ioutil"
	"os"
	"path/filepath"
	"runtime"
	"sort"
	"strings"
)

// This command generate necessary strain information.
type cmdInit struct {
	complete  *bool
	cmdConfig // embed cmdConfig
}

func (cmd *cmdInit) Flags(fs *flag.FlagSet) *flag.FlagSet {
	cmd.cmdConfig.Flags(fs)
	cmd.complete = fs.Bool("complete", false, "only generate completed species map?")
	return fs
}

// Run command.
// It requires
// 1. workspace -- also as output folder;
// 2. reference genome database folder;
// 3. genome reports folder;
// 4. taxonomy database folder.
func (cmd *cmdInit) Run(args []string) {
	// Parse config and settings.
	// Here we will use
	// 1. refBase
	// 2. repBase
	// 3. taxBase
	cmd.ParseConfig()

	// If there is no file containing strain informations,
	// or the one is older than summary.txt,
	// create a new one.
	var strains []strain.Strain
	strains = getStrainInfors(cmd.refBase, cmd.repBase, cmd.taxBase, *cmd.complete)
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

	// Create a specie: []strain map.
	speciesMap := make(map[string][]strain.Strain)
	for _, s := range strains {
		if s.Species != "" {
			if *cmd.complete {
				if strings.Contains(strings.ToLower(s.Status), "complete") {
					speciesMap[s.Species] = append(speciesMap[s.Species], s)
				}
			} else {
				speciesMap[s.Species] = append(speciesMap[s.Species], s)
			}
		}
	}

	speciesNames := []string{}
	for name, _ := range speciesMap {
		speciesNames = append(speciesNames, name)
	}

	sort.Strings(speciesNames)
	yamlFilePath := filepath.Join(*cmd.workspace, "reference_species.yaml")
	w, err = os.Create(yamlFilePath)
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
// after that of prokaryotes.txt in GENOME_REPORTS folder.
func isReferenceStrainsExists(workspace, repBase string) (isExist bool) {
	filePath := filepath.Join(workspace, "reference_strains.json")
	if fi1, err := os.Stat(filePath); err != nil {
		if os.IsNotExist(err) {
			isExist = false
		} else {
			ERROR.Fatalln(err)
		}
	} else {
		prokaryotesFile := filepath.Join(repBase, "prokaryotes.txt")
		fi2, err := os.Stat(prokaryotesFile)
		if err != nil {
			if os.IsNotExist(err) {
				ERROR.Fatalln("Cannot find prokaryotes.txt in genome report directory!")
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

// get strain informations
// from GENOME_REPORTS
func getStrainInfors(refBase, repBase, taxBase string, completed bool) (strains []strain.Strain) {
	// Read prokaryotes strains.
	fileName := "prokaryotes.txt"
	filePath := filepath.Join(repBase, fileName)
	f, err := os.Open(filePath)
	if err != nil {
		ERROR.Panicln(err)
	}
	defer f.Close()
	reportStrains := reports.ReadProkaryotes(f)

	// add genetic code information to each strain.
	taxonMap := taxonomy.ReadTaxas(taxBase)

	jobs := make(chan reports.Strain)
	go func() {
		defer close(jobs)
		for i := 0; i < len(reportStrains); i++ {
			jobs <- reportStrains[i]
		}
	}()

	ncpu := runtime.GOMAXPROCS(0)

	results := make(chan strain.Strain)
	done := make(chan bool)
	for i := 0; i < ncpu; i++ {
		go func() {
			for s := range jobs {
				if s.Path == "-" {
					continue
				}
				s1 := strain.Strain{}
				s1.Name = s.Name
				s1.Path = s.Path
				s1.ProjectId = s.ProjectId
				s1.TaxId = s.TaxId
				s1.Status = s.Status

				taxon, found := taxonMap[s1.TaxId]
				if found {
					s1.GeneticCode = taxon.GeneticCode.Id
					s1.Species = getSpeciesName(s1.TaxId, taxonMap)
				}

				if len(s.Genomes) > 0 {
					for j := 0; j < len(s.Genomes); j++ {
						g1 := s.Genomes[j]
						g2 := genome.Genome{}
						g2.Accession = g1.Accession
						g2.Replicon = "chromosome"
						s1.Genomes = append(s1.Genomes, g2)
					}
				} else {
					if !completed {
						path := filepath.Join(refBase, s.Path)
						genomes := listScaffoldFiles(path)
						for _, g := range genomes {
							g2 := genome.Genome{}
							g2.Accession = g
							g2.Replicon = "chromosome"
							s1.Genomes = append(s1.Genomes, g2)
						}
					}
				}

				results <- s1
			}

			done <- true
		}()
	}

	go func() {
		defer close(results)
		for i := 0; i < ncpu; i++ {
			<-done
		}
	}()

	for s := range results {
		if len(s.Genomes) > 0 && s.Species != "" {
			strains = append(strains, s)
		}
	}

	return
}

func listScaffoldFiles(path string) (genomes []string) {
	fileInfors, err := ioutil.ReadDir(path)
	if err != nil {
		WARN.Println(err)
		return
	}

	scaffoldNameSet := make(map[string]bool)
	for _, fi := range fileInfors {
		if strings.Contains(fi.Name(), "fna.tgz") {
			name := strings.Split(fi.Name(), ".")[0]
			scaffoldNameSet[name] = true
		}
	}

	for name, _ := range scaffoldNameSet {
		genomes = append(genomes, name)
	}

	return
}

// find species for a strain given its taxid.
func getSpeciesName(id string, taxonMap map[string]taxonomy.Taxon) (name string) {
	t := taxonMap[id]

	for {
		if t.Rank == "species" {
			name = t.Name
			break
		} else {
			t = taxonMap[t.Parent]
		}
	}

	return
}
