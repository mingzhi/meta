package main

import (
	"archive/tar"
	"bufio"
	"bytes"
	"compress/gzip"
	"fmt"
	"github.com/mingzhi/biogo/seq"
	"io"
	"os"
	"path/filepath"
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

	for _, strains := range cmd.speciesMap {
		for _, s := range strains {
			if !strings.Contains(strings.ToLower(s.Status), "complete") {
				path := filepath.Join(cmd.refBase, s.Path)
				for _, genome := range s.Genomes {
					acc := genome.Accession
					posMap := mergeFna(acc, path)
					mergeFaa(acc, path)
					mergePtt(acc, path, posMap)
					INFO.Println(path)
				}
			}

		}
	}
}

func mergeFna(genome, path string) map[string]int {
	fileName := genome + ".scaffold.fna.tgz"
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
