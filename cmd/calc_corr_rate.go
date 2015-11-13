package main

import (
	"bufio"
	"bytes"
	"flag"
	"fmt"
	"github.com/mingzhi/biogo/seq"
	"github.com/mingzhi/gomath/stat/correlation"
	"github.com/mingzhi/gomath/stat/desc"
	"github.com/mingzhi/gomath/stat/desc/meanvar"
	"github.com/mingzhi/ncbiftp/seqrecord"
	"github.com/mingzhi/ncbiftp/taxonomy"
	"io"
	"math"
	"os"
	"regexp"
	"strconv"
	"strings"
)

var (
	maxl int
	codonTableId string
	genomeFile string
	pttFile string
	pileFile string
)

func init() {
	flag.IntVar(&maxl, "maxl", 1000, "max length of correlation")
	flag.StringVar(&codonTableId, "codon", "11", "codon table id")
	flag.Parse()
	
	genomeFile = flag.Arg(0)
	pttFile = flag.Arg(1)
	pileFile = flag.Arg(2)
}

func main() {
	s := readGenomeSequence("NC_017501" + ".fna")[0]

	ptts := readPtt("NC_017501" + ".ptt")

	// determine codon table.
	codonTableMap := taxonomy.GeneticCodes()
	codonTable := codonTableMap[codonTableId]
	profile := profileGenome(s.Seq, ptts, gc)

	filename := "mapped.pileup"
	f, err := os.Open(filename)
	if err != nil {
		panic(err)
	}
	defer f.Close()

	reader := bufio.NewReader(f)
	pis := decode(reader)
	totalPis := divivePis(pis, 100)
	totalCorrs := [][]float64{}
	for i := 0; i < len(totalPis); i++ {
		corrs, _ := calcCov(totalPis[i], profile, ThirdPos)
		totalCorrs = append(totalCorrs, corrs)
	}
	corrs, vars, ns := collect(totalCorrs)

	w, err := os.Create("test.txt")
	if err != nil {
		panic(err)
	}
	defer w.Close()
	for i := 0; i < len(corrs); i++ {
		w.WriteString(fmt.Sprintf("%d\t%g\t%g\t%d\n", i, corrs[i], vars[i], ns[i]))
	}
}

func profileGenome(genomeFile, pttFile string) (profile []byte) {
	profile = make([]byte, len(s))
	for i := 0; i < len(profile); i++ {
		profile[i] = NonCoding
	}

	for _, ptt := range ptts {
		// Prepare nucleotide sequence,
		// we need it for determine 4-fold codons.
		var nucl []byte
		if ptt.Loc.To >= ptt.Loc.From {
			nucl = s[ptt.Loc.From-1 : ptt.Loc.To]
		} else {
			// skip genes across boundary.
			continue
		}

		if ptt.Loc.Strand == "-" {
			nucl = seq.Complement(seq.Reverse(nucl))
		}

		prof := make([]byte, len(nucl))
		for j, _ := range nucl {
			switch (j + 1) % 3 {
			case 1:
				prof[j] = FirstPos
			case 2:
				prof[j] = SecondPos
			case 0:
				codon := nucl[j-2 : j+1]
				if gc.FFCodons[string(codon)] {
					prof[j] = FourFold
				} else {
					prof[j] = ThirdPos
				}
			}
		}

		if ptt.Loc.Strand == "-" {
			prof = seq.Reverse(prof)
		}

		for j, p := range prof {
			index := ptt.Loc.From - 1 + j
			if profile[index] == NonCoding {
				profile[index] = p
			} else {
				profile[index] = Undefined
			}
		}
	}
	return
}

type SNP struct {
	Genome    string
	Position  int
	RefBase   byte
	Number    int
	ReadBases []byte
	BaseQuals []int
}

type Pi struct {
	Position int
	Pi       float64
}

func decode(reader *bufio.Reader) []Pi {
	total := 0
	n := 0
	ks := desc.NewMean()
	pis := []Pi{}
	for {
		line, err := reader.ReadString('\n')
		if err != nil {
			if err != io.EOF {
				panic(err)
			} else {
				break
			}
		}
		trimedLine := strings.TrimSpace(line)
		terms := strings.Split(trimedLine, "\t")
		snp := SNP{}
		snp.Genome = terms[0]
		snp.Position = str2int(terms[1])
		snp.RefBase = terms[2][0]
		snp.Number = str2int(terms[3])
		snp.ReadBases = decodeReadBases(terms[4])
		snp.BaseQuals = decodeBaseQual(terms[5])
		if snp.Number != len(snp.ReadBases) || snp.Number != len(snp.BaseQuals) {
			fmt.Printf("%d\t%s\n", snp.Number, string(snp.ReadBases))
		}
		snp.ReadBases = filterBases(snp.ReadBases, snp.BaseQuals, 30)
		if snp.Number >= 10 {
			pi := calcPi(snp.ReadBases)
			if !math.IsNaN(pi) {
				n++
				ks.Increment(pi)
				pis = append(pis, Pi{Position: snp.Position, Pi: pi})
			}
		}
		total++
	}
	fmt.Printf("%d\t%d\n", n, total)
	fmt.Printf("%f\t%d\n", ks.GetResult(), ks.GetN())

	return pis
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

func calcPi(s []byte) (pi float64) {
	us := bytes.ToUpper(s)
	total := 0
	n := 0

	numBase := 0
	for i := 0; i < len(us); i++ {
		if us[i] != '*' {
			numBase++
		}
	}
	if numBase <= 10 {
		pi = math.NaN()
		return
	}

	for i := 0; i < len(us); i++ {
		for j := i + 1; j < len(us); j++ {
			if us[i] != '*' && us[j] != '*' {
				if us[i] != us[j] {
					n++
				}
				total++
			}
		}
	}

	pi = float64(n) / float64(total)

	return
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

func calcCov(pis []Pi, profile []byte, pos byte) (covs []float64, n []int) {
	maxl := 300
	corrs := make([]*correlation.BivariateCovariance, maxl)
	for i := 0; i < maxl; i++ {
		corrs[i] = correlation.NewBivariateCovariance(false)
	}

	for i := 0; i < len(pis); i++ {
		pos1 := profile[pis[i].Position]
		if pos == ThirdPos {
			if pos1 == FourFold {
				pos1 = ThirdPos
			}
		}
		if pos1 == pos {
			for j := i; j < len(pis); j++ {
				pos2 := profile[pis[j].Position]
				if pos == ThirdPos {
					if pos2 == FourFold {
						pos2 = ThirdPos
					}
				}
				if pos2 == pos {
					l := pis[j].Position - pis[i].Position
					if l < maxl {
						corrs[l].Increment(pis[i].Pi, pis[j].Pi)
					} else {
						break
					}
				}
			}
		}

	}

	for i := 0; i < maxl; i++ {
		covs = append(covs, corrs[i].GetResult())
		n = append(n, corrs[i].GetN())
	}

	return
}

func collect(corrs [][]float64) (means []float64, vars []float64, ns []int) {
	mvs := make([]*meanvar.MeanVar, len(corrs[0]))
	for i := 0; i < len(mvs); i++ {
		mvs[i] = meanvar.New()
	}

	for i := 0; i < len(mvs); i++ {
		for j := 0; j < len(corrs); j++ {
			v := corrs[j][i]
			if !math.IsNaN(v) {
				mvs[i].Increment(v)
			}
		}
	}

	means = make([]float64, len(mvs))
	vars = make([]float64, len(mvs))
	ns = make([]int, len(mvs))
	for i := 0; i < len(mvs); i++ {
		means[i] = mvs[i].Mean.GetResult()
		ns[i] = mvs[i].Mean.GetN()
		vars[i] = mvs[i].Var.GetResult()
	}

	return
}

func divivePis(pis []Pi, num int) [][]Pi {
	size := len(pis) / num
	totalPis := [][]Pi{}
	to := size
	chose := []Pi{}
	for i := 0; i < len(pis); i++ {
		if i >= to {
			to += size
			totalPis = append(totalPis, chose)
			chose = []Pi{}
		}
		chose = append(chose, pis[i])
	}
	totalPis = append(totalPis, chose)
	return totalPis
}

// Read genome sequence using biogo.
func readGenomeSequence(filename string) []*seq.Sequence {
	// open file stream.
	f, err := os.Open(filename)
	if err != nil {
		panic(err)
	}
	defer f.Close()

	rd := seq.NewFastaReader(f)
	sequences, err := rd.ReadAll()
	if err != nil {
		panic(err)
	}

	return sequences
}

// Read ptt file.
func readPtt(filename string) []seqrecord.Ptt {
	pttReader := seqrecord.NewPttFile(filename)
	ptts := pttReader.ReadAll()

	return ptts
}

const (
	NonCoding byte = '0'
	FirstPos  byte = '1'
	SecondPos byte = '2'
	ThirdPos  byte = '3'
	FourFold  byte = '4'
	Undefined byte = '5'
)
