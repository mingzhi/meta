package main

import (
	"archive/tar"
	"bufio"
	"bytes"
	"compress/gzip"
	"fmt"
	"github.com/mingzhi/biogo/seq"
	"github.com/mingzhi/meta/strain"
	"io"
	"os"
	"path/filepath"
	"runtime"
	"strconv"
	"strings"
)

// This command merge scaffolds.
type cmdScaffoldMerge struct {
	cmdConfig
}

func (cmd *cmdScaffoldMerge) Run(args []string) {
	// Parse config and settings.
	cmd.ParseConfig()
	// Load species strain maps.
	cmd.LoadSpeciesMap()
	completed := []string{
		"Complete",
		"Complete Genome",
		"Chromosome with gaps",
		"Gapless Chromosome",
	}

	type job struct {
		acc, path string
	}
	jobs := make(chan job)
	go func() {
		for _, strains := range cmd.speciesMap {
			for _, s := range strains {
				if !stringInSlice(strings.TrimSpace(s.Status), completed) {
					path := filepath.Join(cmd.refBase, s.Path)
					for _, genome := range s.Genomes {
						acc := genome.RefAcc()
						jobs <- job{acc, path}
					}
				}
			}
		}
		close(jobs)
	}()

	done := make(chan bool)
	ncpu := runtime.GOMAXPROCS(0)
	for i := 0; i < ncpu; i++ {
		go func() {
			for job := range jobs {
				path := job.path
				acc := job.acc
				posMap := mergeFna(acc, path)
				if posMap != nil {
					mergeFaa(acc, path)
					mergePtt(acc, path, posMap)
				}
			}
			done <- true
		}()
	}

	for i := 0; i < ncpu; i++ {
		<-done
	}
}

// Load species map.
func (cmd *cmdScaffoldMerge) LoadSpeciesMap() {
	// If reference_strains exists, read it,
	// else generate it.
	var strains []strain.Strain
	if isReferenceStrainsExists(*cmd.workspace, cmd.repBase) {
		strains = cmd.ReadReferenceStrains()
	} else {
		ERROR.Fatalln("Can not find reference_strains.json, please run meta init first!")
	}

	// Make a strain map,
	// strain.path: strain
	strainMap := make(map[string]strain.Strain)
	for _, strain := range strains {
		strainMap[strain.Path] = strain
	}

	// Read species input yaml file.
	inputMap := cmd.ReadSpeciesFile()

	// Load species map.
	cmd.speciesMap = make(map[string][]strain.Strain)
	for prefix, strainPaths := range inputMap {
		for _, strainPath := range strainPaths {
			strain, found := strainMap[strainPath]
			if found {
				cmd.speciesMap[prefix] = append(cmd.speciesMap[prefix], strain)
			} else {
				WARN.Printf("Cannot find %s\n", strainPath)
			}
		}
	}
}

func stringInSlice(a string, list []string) bool {
	for _, b := range list {
		if b == a {
			return true
		}
	}
	return false
}

func mergeFna(genome, path string) map[string]int {
	fileName := genome + ".scaffold.fna.tgz"
	filePath := filepath.Join(path, fileName)
	f, err := os.Open(filePath)
	if err != nil {
		WARN.Println(err)
		return nil
	}
	defer f.Close()

	// handle gzipped.
	rd, err := gzip.NewReader(f)
	if err != nil {
		ERROR.Panic(err)
	}
	defer rd.Close()

	tarReader := tar.NewReader(rd)

	masks := bytes.Repeat([]byte{'N'}, 1000)

	gs := seq.Sequence{}
	posMap := make(map[string]int)
	// extract files.
	for {
		header, err := tarReader.Next()
		if err != nil {
			if err == io.EOF {
				break
			} else {
				ERROR.Panicln(err)
			}
		}

		acc := strings.Split(header.Name, ".")[0]
		posMap[acc] = len(gs.Seq)

		// read fasta sequence.
		fastaReader := seq.NewFastaReader(tarReader)
		sequences, err := fastaReader.ReadAll()
		if err != nil {
			ERROR.Panicln(err)
		}
		for _, s := range sequences {
			gs.Seq = append(gs.Seq, s.Seq...)
		}

		gs.Seq = append(gs.Seq, masks...)
	}

	outName := genome + ".fna"
	outPath := filepath.Join(path, outName)
	out, err := os.Create(outPath)
	if err != nil {
		ERROR.Panicln(err)
	}
	defer out.Close()

	out.WriteString(fmt.Sprintf(">ref|%s\n", genome))
	out.Write(gs.Seq)

	return posMap
}

func mergeFaa(genome, path string) {
	fileName := genome + ".scaffold.faa.tgz"
	filePath := filepath.Join(path, fileName)
	f, err := os.Open(filePath)
	if err != nil {
		ERROR.Panic(err)
	}
	defer f.Close()

	// handle gzipped.
	rd, err := gzip.NewReader(f)
	if err != nil {
		ERROR.Panic(err)
	}
	defer rd.Close()

	tarReader := tar.NewReader(rd)

	outName := genome + ".faa"
	outPath := filepath.Join(path, outName)
	out, err := os.Create(outPath)
	if err != nil {
		ERROR.Panicln(err)
	}
	defer out.Close()

	// extract files.
	for {
		_, err := tarReader.Next()
		if err != nil {
			if err == io.EOF {
				break
			} else {
				ERROR.Panicln(err)
			}
		}
		io.Copy(out, tarReader)
	}
}

func mergePtt(genome, path string, posMap map[string]int) {
	fileName := genome + ".scaffold.ptt.tgz"
	filePath := filepath.Join(path, fileName)
	f, err := os.Open(filePath)
	if err != nil {
		ERROR.Panic(err)
	}
	defer f.Close()

	// handle gzipped.
	rd, err := gzip.NewReader(f)
	if err != nil {
		ERROR.Panic(err)
	}
	defer rd.Close()

	tarReader := tar.NewReader(rd)

	outName := genome + ".ptt"
	outPath := filepath.Join(path, outName)
	out, err := os.Create(outPath)
	if err != nil {
		ERROR.Panicln(err)
	}
	defer out.Close()

	// extract files.
	for {
		header, err := tarReader.Next()
		if err != nil {
			if err == io.EOF {
				break
			} else {
				ERROR.Panicln(err)
			}
		}
		acc := strings.Split(header.Name, ".")[0]
		pos, found := posMap[acc]
		if found {
			bufReader := bufio.NewReader(tarReader)
			for {
				line, err := bufReader.ReadString('\n')
				if err != nil {
					if err == io.EOF {
						break
					} else {
						ERROR.Panicln(err)
					}
				}

				terms := strings.Split(line, "\t")
				if len(terms) > 5 && strings.Contains(terms[0], "..") {
					startend := strings.Split(terms[0], "..")
					start, _ := strconv.Atoi(startend[0])
					end, _ := strconv.Atoi(startend[1])
					fields := []string{fmt.Sprintf("%d..%d", start+pos, end+pos)}
					fields = append(fields, terms[1:]...)
					out.WriteString(strings.Join(fields, "\t"))
				}
			}
		}
	}
}
