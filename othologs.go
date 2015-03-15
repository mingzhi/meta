package meta

import (
	"github.com/mingzhi/ncbiutils"
	"path/filepath"
	"runtime"
)

var GCMAP map[string]*ncbiutils.GeneticCode

// Return sequence records for each ortholog clusters.
func FindOrthologs(strains []Strain, dir string, clusters [][]string) []ncbiutils.SeqRecords {
	// Load sequence records for each genome.
	ncpu := runtime.GOMAXPROCS(0)
	jobs := make(chan []string)
	go func() {
		for i := 0; i < len(strains); i++ {
			s := strains[i]
			d := filepath.Join(dir, s.Path)
			for _, g := range s.Genomes {
				jobs <- []string{g.Accession, d, s.GeneticCode}
			}
		}
		close(jobs)
	}()

	done := make(chan bool)
	results := make(chan ncbiutils.SeqRecords)
	for i := 0; i < ncpu; i++ {
		go func() {
			for job := range jobs {
				g, d, gc := job[0], job[1], job[2]
				records := ncbiutils.ReadSeqRecords(g, d, gc)
				results <- records
			}
			done <- true
		}()
	}

	go func() {
		for i := 0; i < ncpu; i++ {
			<-done
		}
		close(results)
	}()

	recMap := make(map[string]ncbiutils.SeqRecord)
	for rec := range results {
		for _, r := range rec {
			recMap[r.Id+"|"+r.Genome] = r
		}
	}

	oGroups := []ncbiutils.SeqRecords{}
	for _, cluster := range clusters {
		grp := ncbiutils.SeqRecords{}
		for _, s := range cluster {
			r, found := recMap[s]
			if found {
				grp = append(grp, r)
			}
		}
		oGroups = append(oGroups, grp)
	}

	return oGroups
}
