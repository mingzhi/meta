package main

import (
	"bufio"
	"bytes"
	"encoding/json"
	"flag"
	"io"
	"log"
	"math"
	"os"
	"regexp"
	"strconv"
	"strings"
)

func main() {
	var filename string
	var outfile string
	flag.Parse()
	if flag.NArg() < 1 {
		log.Fatalln("Usage: go run pileup2pi.go <pileup file> <output file>")
	}
	filename = flag.Arg(0)
	outfile = flag.Arg(1)

	f, err := os.Open(filename)
	if err != nil {
		log.Fatalln(err)
	}
	defer f.Close()

	// Create output file.
	w, err := os.Create(outfile)
	if err != nil {
		log.Fatalln(err)
	}
	defer w.Close()
	encoder := json.NewEncoder(w)

	snpChan := DecodePileup(f)
	for s := range snpChan {
		pi := CalcPi(s)
		if !math.IsNaN(pi.Pi) {
			err := encoder.Encode(pi)
			if err != nil {
				log.Fatalln(err)
			}
		}
	}
}

type Pi struct {
	Genome   string
	Position int
	Pi       float64
	Num      int
	Alleles  map[string]int
}

// Calc Pi from each snp.
func CalcPi(snp SNP) (pi Pi) {
	pi.Genome = snp.Genome
	pi.Position = snp.Position
	pi.Pi, pi.Num, pi.Alleles = calcPi(snp.ReadBases)
	return pi
}

func calcPi(bases []byte) (pi float64, n int, m map[string]int) {
	// convert bases to upper case.
	upperBases := bytes.ToUpper(bases)

	m = make(map[string]int)
	for i := 0; i < len(upperBases); i++ {
		if upperBases[i] != '*' {
			b := string(upperBases[i])
			m[b]++
		}
	}

	total := 0
	nums := []int{}
	for _, n := range m {
		total += n
		nums = append(nums, n)
	}

	cross := 0
	for i := 0; i < len(nums); i++ {
		for j := i + 1; j < len(nums); j++ {
			cross += nums[i] * nums[j]
		}
	}

	pi = float64(cross) / float64(total*(total-1)/2)
	n = total
	return
}

// SNP contains single nucleotide polymorphisms.
type SNP struct {
	Genome    string
	Position  int // 1-coordinate system [1 - N]
	RefBase   byte
	Number    int
	ReadBases []byte
	BaseQuals []int
}

func DecodePileup(r io.Reader) (snpChan chan SNP) {
	snpChan = make(chan SNP)
	go func() {
		defer close(snpChan)
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
			if snp.Number > 0 {
				snp.ReadBases = decodeReadBases(terms[4], snp.RefBase)
				snp.BaseQuals = decodeBaseQual(terms[5])
				if snp.Number != len(snp.ReadBases) || snp.Number != len(snp.BaseQuals) {
					log.Printf("SNP number and length of bases did not matched: Position %d, Number %d, Bases %s\n", snp.Position, snp.Number, string(snp.ReadBases))
				} else {
					snp.ReadBases = filterBases(snp.ReadBases, snp.BaseQuals, 30)
					snpChan <- snp
				}
			}
		}
	}()

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

func decodeReadBases(s string, ref byte) []byte {
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
		log.Printf("Bases still contains number: %s\n", s)
	}

	bases := []byte{}
	for i := 0; i < len(s); i++ {
		b := s[i]
		if b == '.' || b == ',' {
			bases = append(bases, ref)
		} else {
			bases = append(bases, b)
		}
	}

	bases = bytes.ToUpper(bases)

	return bases
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
