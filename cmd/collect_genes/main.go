package main

import (
	"bufio"
	"encoding/json"
	"fmt"
	"io"
	"log"
	"os"
	"sort"
	"strings"

	"github.com/alecthomas/kingpin"
)

func main() {
	var corrFile string
	var outfile string
	var geneFile string
	var sampleFile string
	var byGene bool
	app := kingpin.New("collect_genes", "Calculate correlation across multiple samples")
	app.Version("v0.1")

	outFileArg := app.Arg("out-file", "output file").Required().String()
	corrFileArg := app.Flag("corr-res-file", "corr results file").Default("").String()
	geneFileFlag := app.Flag("gene-file", "gene file").Default("").String()
	sampleFileFlag := app.Flag("sample-file", "sample file").Default("").String()
	byGeneFlag := app.Flag("by-gene", "by gene").Default("false").Bool()
	kingpin.MustParse(app.Parse(os.Args[1:]))
	corrFile = *corrFileArg
	outfile = *outFileArg
	geneFile = *geneFileFlag
	sampleFile = *sampleFileFlag
	byGene = *byGeneFlag

	var geneSet map[string]bool
	if geneFile != "" {
		geneSet = make(map[string]bool)
		lines := readLines(geneFile)
		for _, line := range lines {
			gene := strings.Split(line, "\t")[0]
			geneSet[gene] = true
		}
	}
	var samples []string
	if sampleFile != "" {
		lines := readLines(sampleFile)
		for _, line := range lines {
			samples = append(samples, strings.TrimSpace(line))
		}
	} else {
		samples = append(samples, corrFile)
	}

	collectorMap := make(map[string]*Collector)
	if !byGene {
		collectorMap["all"] = NewCollector()
	}
	for _, sampleFile := range samples {
		corrChan := readCorrResults(sampleFile)
		for corrResults := range corrChan {
			geneID := corrResults.GeneID
			if geneFile != "" {
				if !geneSet[geneID] {
					continue
				}
			}

			collectorMap["all"].Add(corrResults)
		}
	}

	w, err := os.Create(outfile)
	if err != nil {
		log.Panic(err)
	}
	defer w.Close()

	// sort by gene id
	var geneIDs []string
	for geneID := range collectorMap {
		geneIDs = append(geneIDs, geneID)
	}
	sort.Strings(geneIDs)

	w.WriteString("l,m,v,n,t,g\n")
	for _, geneID := range geneIDs {
		collector := collectorMap[geneID]
		results := collector.Results()
		for _, res := range results {
			w.WriteString(fmt.Sprintf("%d,%g,%g,%d,%s,%s\n",
				res.Lag, res.Value, res.Variance, res.Count, res.Type, geneID))
		}
	}
}

func readSamples(filename string) []string {
	f, err := os.Open(filename)
	if err != nil {
		log.Panic(err)
	}
	defer f.Close()

	rd := bufio.NewReader(f)
	var results []string
	for {
		line, err := rd.ReadString('\n')
		if err != nil {
			if err != io.EOF {
				log.Panic(err)
			}
			break
		}
		results = append(results, strings.TrimSpace(line))
	}
	return results
}

func readCorrResults(filename string) chan CorrResults {
	c := make(chan CorrResults)
	go func() {
		defer close(c)
		f, err := os.Open(filename)
		if err != nil {
			log.Panic(err)
		}
		defer f.Close()
		decoder := json.NewDecoder(f)
		for {
			var rec CorrResults
			if err := decoder.Decode(&rec); err != nil {
				if err != io.EOF {
					log.Panic(err)
				}
				break
			}
			c <- rec
		}
	}()
	return c
}

func checkFiles(samples []string, appendix string) []string {
	var results []string
	for _, sample := range samples {
		filename := sample + appendix
		if _, err := os.Stat(filename); !os.IsNotExist(err) {
			results = append(results, sample)
		}
	}
	return results
}

// readLines return all trimmed lines.
func readLines(filename string) []string {
	f, err := os.Open(filename)
	if err != nil {
		log.Panic(err)
	}
	defer f.Close()

	rd := bufio.NewReader(f)
	var lines []string
	for {
		line, err := rd.ReadString('\n')
		if err != nil {
			if err != io.EOF {
				log.Panic(err)
			}
			break
		}
		lines = append(lines, strings.TrimSpace(line))
	}
	return lines
}
