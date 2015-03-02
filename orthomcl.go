package meta

// Perform ortholog clustering using reciprocal best hit and MCL.

import (
	"bufio"
	"bytes"
	"fmt"
	"github.com/mingzhi/hgt/taxa"
	"io"
	"io/ioutil"
	"log"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"strings"
)

// OrthoMCL identifies ortholog groups for a group of closed related strains.
// First, it uses usearch to find the best reciprocal top hits
// for each pair of genomes.
// Then it uses MCL to obtain ortholog clusters.
// strains: a array of strains.
// dir: reference genome folder, containing their genome sequences.
func OrthoMCl(strains []Strain, dir string) (clusters [][]string) {
	// Check and prepare blast database.
	// Remove those strains that do not have protein sequences.
	selectStrains := []Strain{}
	for _, s := range strains {
		for _, g := range s.Genomes {
			acc := cleanAccession(g.Accession)
			f := filepath.Join(dir, s.Path, acc+".faa")
			if _, err := os.Stat(f); err == nil {
				// index the sequence if have not.
				if !IsUsearchDBExist(f) {
					UsearchMakeUDB(f)
				}
				selectStrains = append(selectStrains, s)
			} else if os.IsNotExist(err) {
				Warn.Printf("%s is not exist!\n", f)
			}
		}
	}

	// Perform all against all usearch.
	pairs := AllAgainstAll(selectStrains, dir)

	// Get ortholog clusters using MCL
	clusters = MCL(pairs)

	Info.Printf("Total clusters: %d\n", len(clusters))
	n := 0.0
	for _, cls := range clusters {
		n += float64(len(cls))
	}
	Info.Printf("Average cluster size: %.1f\n", n/float64(len(clusters)))

	return
}

type Pair struct {
	A, B  string
	Score float64
}

// AllAgainstAll perform all against all usearch
// for the reciprocal top hits for each pair of genomes.
func AllAgainstAll(strains []Strain, dir string) []Pair {
	ncpu := runtime.GOMAXPROCS(0)
	// Prepare jobs.
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

	done := make(chan bool)      // signal channel.
	results := make(chan []Pair) // result channel.
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
			done <- true
		}()
	}

	go func() {
		for i := 0; i < ncpu; i++ {
			<-done
		}
		close(results)
	}()

	allPairs := []Pair{}
	n := 0
	for pairs := range results {
		allPairs = append(allPairs, pairs...)
		n++
	}
	Info.Printf("Total orthologous pairs: %d\n", len(allPairs))
	Info.Printf("Total pairs of genomes: %d, average %.1f pairs of orthologs\n",
		n, float64(len(allPairs))/float64(n))
	return allPairs
}

// Find reciprocal top hits for a pair of genomes.
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

func MCL(pairs []Pair) (clusters [][]string) {
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

	cmd := exec.Command("mcl", "-", "--abc", "-o", f.Name())
	cmd.Stdin = content

	stderr := new(bytes.Buffer)
	cmd.Stderr = stderr

	if err := cmd.Run(); err != nil {
		panic(string(stderr.Bytes()))
	}

	clusters = readMCL(f.Name())
	return
}

func readMCL(filename string) (clusters [][]string) {
	f, err := os.Open(filename)
	if err != nil {
		panic(err)
	}
	defer f.Close()
	r := bufio.NewReader(f)
	for {
		l, err := r.ReadString('\n')
		if err != nil {
			if err != io.EOF {
				panic(err)
			}
			break
		}
		l = strings.TrimSpace(l)
		genes := strings.Split(l, "\t")
		clusters = append(clusters, genes)
	}

	return
}

func LoadSeqRecords(strains []Strain, dir string, gcMap map[string]*taxa.GeneticCode) {

}

func cleanAccession(name string) string {
	return strings.Split(name, ".")[0]
}
