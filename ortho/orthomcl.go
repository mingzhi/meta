package ortho

// Perform ortholog clustering using reciprocal best hit and MCL.

import (
	"bufio"
	"bytes"
	"fmt"
	"github.com/mingzhi/meta/strain"
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
//
// strains: a array of strains.
// dir: reference genome folder, containing their genome sequences.
func OrthoMCl(strains []strain.Strain, dir string) (clusters [][]string) {
	// Check and prepare blast database.
	// Remove those strains that do not have protein sequences.
	selectStrains := []strain.Strain{}
	for _, s := range strains {
		for _, g := range s.Genomes {
			f := filepath.Join(dir, s.Path, g.RefAcc()+".faa")
			if _, err := os.Stat(f); err == nil {
				// index the sequence if have not.
				if !IsUsearchDBExist(f) {
					UsearchMakeUDB(f)
				}
				selectStrains = append(selectStrains, s)
			} else {
				panic(err)
			}
		}
	}

	// Perform all against all usearch.
	pairs := AllAgainstAll(selectStrains, dir)

	// Get ortholog clusters using MCL
	clusters = MCL(pairs)

	return
}

type Pair struct {
	A, B  string
	Score float64
}

// AllAgainstAll perform all against all usearch
// for the reciprocal top hits for each pair of genomes.
func AllAgainstAll(strains []strain.Strain, dir string) []Pair {
	ncpu := runtime.GOMAXPROCS(0)
	// Prepare jobs.
	jobs := make(chan []strain.Strain, ncpu)
	go func() {
		for i := 0; i < len(strains); i++ {
			a := strains[i]
			for j := i + 1; j < len(strains); j++ {
				b := strains[j]
				pair := []strain.Strain{a, b}
				jobs <- pair

			}
		}
		close(jobs)
	}()

	done := make(chan bool)    // signal channel.
	results := make(chan Pair) // result channel.
	for i := 0; i < ncpu; i++ {
		go func() {
			for pair := range jobs {
				a := pair[0]
				b := pair[1]
				for _, genomeA := range a.Genomes {
					accA := genomeA.RefAcc()
					for _, genomeB := range b.Genomes {
						accB := genomeB.RefAcc()
						fileNameA := filepath.Join(dir, a.Path, accA+".faa")
						fileNameB := filepath.Join(dir, b.Path, accB+".faa")
						pairs := ReciprocalBestHits(fileNameA, fileNameB)
						for _, p := range pairs {
							newP := Pair{
								A: p.A + "|" + accA,
								B: p.B + "|" + accB,
							}
							results <- newP
						}
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

	for pair := range results {
		allPairs = append(allPairs, pair)
	}

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
