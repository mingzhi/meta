package main

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"regexp"
	"strconv"
	"strings"
)

func DecodePileup(r io.Reader, snpChan chan SNP) {
	reader := bufio.NewReader(r)
	for {
		line, endOfFile := readLine(reader)
		if endOfFile {
			break
		}

		// tab separated line.
		trimedLine := strings.TrimSpace(line)
		terms := strings.Split(trimedLine, "\t")

		// create a SNP
		snp := SNP{}
		snp.Genome = terms[0]
		snp.Position = str2int(terms[1])
		snp.RefBase = terms[2][0]
		snp.Number = str2int(terms[3])
		snp.ReadBases = decodeReadBases(terms[4])
		snp.BaseQuals = decodeBaseQual(terms[5])
		if snp.Number != len(snp.ReadBases) || snp.Number != len(snp.BaseQuals) {
			fmt.Printf("%d\t%s\n", snp.Number, string(snp.ReadBases))
		} else {
			snp.ReadBases = filterBases(snp.ReadBases, snp.BaseQuals, 30)
			snpChan <- snp
		}
	}

	return
}

func readFile(filename string) *os.File {
	f, err := os.Open(filename)
	if err != nil {
		panic(err)
	}
	return f
}

func readLine(reader *bufio.Reader) (line string, endOfFile bool) {
	l, err := reader.ReadString('\n')
	if err != nil {
		if err != io.EOF {
			panic(err)
		} else {
			endOfFile = true
		}
	}
	line = l

	return
}

func decodeReadBases(s string) []byte {
	r := regexp.MustCompile("\\^.")
	s = r.ReplaceAllString(s, "")
	s = strings.Replace(s, "$", "", -1)

	r2 := regexp.MustCompile("[\\+-][0-9]+")
	if r2.MatchString(s) {
		insertNumbers := r2.FindAllString(s, -1)
		insertIndex := r2.FindAllStringIndex(s, -1)
		deletedPositions := make(map[int]bool)
		for i := 0; i < len(insertNumbers); i++ {
			start := insertIndex[i][0]
			n := str2int(insertNumbers[i][1:])
			for j := start; j < start+2+n; j++ {
				deletedPositions[j] = true
			}
		}
		bs := []byte{}
		for i := 0; i < len(s); i++ {
			if !deletedPositions[i] {
				bs = append(bs, s[i])
			}
		}
		s = string(bs)
	}

	r3 := regexp.MustCompile("[0-9]+")
	if r3.MatchString(s) {
		fmt.Println(s)
	}

	return []byte(s)
}

func decodeBaseQual(s string) []int {
	scores := []int{}
	for i := 0; i < len(s); i++ {
		scores = append(scores, int(s[i]))
	}
	return scores
}

func str2int(s string) int {
	i, err := strconv.Atoi(s)
	if err != nil {
		panic(err)
	}
	return i
}

func filterBases(bases []byte, quals []int, cutoff int) []byte {
	bases1 := []byte{}
	for i := 0; i < len(bases); i++ {
		if quals[i] < cutoff {
			bases1 = append(bases1, '*')
		} else {
			bases1 = append(bases1, bases[i])
		}
	}

	return bases1
}
