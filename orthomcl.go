package meta

import (
	"bytes"
	"fmt"
	"io/ioutil"
	"log"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"strings"
)

func findOrthologs(strains []Strain, dir string) (clusters [][]string) {
	// Check and prepare blast database.
	selectStrains := []Strain{}
	for _, s := range strains {
		for _, g := range s.Genomes {
			f := filepath.Join(dir, s.Path, cleanAccession(g.Accession)+".faa")
			if _, err := os.Stat(f); err == nil {
				if !IsUsearchDBExist(f) {
					UsearchMakeUDB(f)
				}
				selectStrains = append(selectStrains, s)
			} else if os.IsNotExist(err) {
				Warn.Printf("%s is not exist!\n", f)
			}
		}
	}

	// perform all against all blast
	pairs := AllAgainstAll(selectStrains, dir)

	// get ortholog clusters using MCL
	MCL(pairs)

	return
}

type Pair struct {
	A, B string
}

func AllAgainstAll(strains []Strain, dir string) []Pair {
	ncpu := runtime.GOMAXPROCS(0)
	jobs := make(chan []Strain, ncpu)
	go func() {
		for i := 0; i < len(strains); i++ {
			a := strains[i]
			for j := i + 1; j < len(strains); j++ {
				b := strains[j]
				pair := []Strain{a, b}
				jobs <- pair

			}
		}
		close(jobs)
	}()

	done := make(chan bool)
	results := make(chan []Pair)
	for i := 0; i < ncpu; i++ {
		go func() {
			for pair := range jobs {
				a := pair[0]
				b := pair[1]
				for _, genomeA := range a.Genomes {
					accA := cleanAccession(genomeA.Accession)
					for _, genomeB := range b.Genomes {
						accB := cleanAccession(genomeB.Accession)
						fileNameA := filepath.Join(dir, a.Path, accA+".faa")
						fileNameB := filepath.Join(dir, b.Path, accB+".faa")
						pairs := ReciprocalBestHits(fileNameA, fileNameB)
						newPairs := []Pair{}
						for _, p := range pairs {
							newP := Pair{
								A: p.A + "|" + accA,
								B: p.B + "|" + accB,
							}
							newPairs = append(newPairs, newP)
						}
						results <- newPairs
					}
				}
			}
		}()
	}

	go func() {
		for i := 0; i < ncpu; i++ {
			<-done
		}
		close(results)
	}()

	allPairs := []Pair{}
	for pairs := range results {
		allPairs = append(allPairs, pairs...)
	}

	return allPairs
}

func ReciprocalBestHits(a, b string) []Pair {
	hits1 := UsearchGlobal(a, b)
	m1 := hit2Map(hits1)
	hits2 := UsearchGlobal(b, a)
	m2 := hit2Map(hits2)
	pairs := []Pair{}
	for q, s := range m1 {
		s2, found := m2[s]
		if found && q == s2 {
			pairs = append(pairs, Pair{A: q, B: s})
		}
	}

	return pairs
}

func hit2Map(hits []Hit) map[string]string {
	m := make(map[string]string)
	for _, h := range hits {
		_, found := m[h.QSeqid]
		if !found {
			m[h.QSeqid] = h.SSeqid
		}
	}

	return m
}

func MCL(pairs []Pair) {
	content := new(bytes.Buffer)
	for _, p := range pairs {
		content.WriteString(fmt.Sprintf("%s %s 1\n", p.A, p.B))
	}

	f, err := ioutil.TempFile("", "mcl")
	if err != nil {
		log.Panic(err)
	}
	f.Close()
	defer os.Remove(f.Name())

	cmd := exec.Command("mcl", "-", "--abc", "-o", f.Name()+".mcl")
	cmd.Stdin = content

	stderr := new(bytes.Buffer)
	cmd.Stderr = stderr

	if err := cmd.Run(); err != nil {
		panic(string(stderr.Bytes()))
	}

	return
}

func cleanAccession(name string) string {
	return strings.Split(name, ".")[0]
}
