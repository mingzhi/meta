package main

import (
	"flag"
	"fmt"
	"github.com/mingzhi/gomath/stat/desc/meanvar"
	"github.com/mingzhi/ncbiftp/taxonomy"
	"math"
	"os"
	"runtime"
)

var (
	genomeFile   string
	pttFile      string
	pileupFile   string
	outFile      string
	maxl         int
	pos          int
	codonTableId string
)

func init() {
	flag.IntVar(&maxl, "maxl", 1500, "max length of correlation")
	flag.IntVar(&pos, "pos", 3, "codon position")
	flag.StringVar(&codonTableId, "code", "11", "codon table id")
	flag.Parse()
	if flag.NArg() < 4 {
		fmt.Println("Usage: calc_corr_rate <genome file> <ptt file> <pileup file> <out file>")
		os.Exit(1)
	}
	genomeFile = flag.Arg(0)
	pttFile = flag.Arg(1)
	pileupFile = flag.Arg(2)
	outFile = flag.Arg(3)

	runtime.GOMAXPROCS(runtime.NumCPU())
}

func main() {
	codonTable := taxonomy.GeneticCodes()[codonTableId]
	profile := ProfileGenome(genomeFile, pttFile, codonTable)
	fmt.Println("Finish generating profile!")

	f := readFile(pileupFile)
	snpChan := make(chan SNP)
	defer f.Close()
	go func() {
		defer close(snpChan)
		DecodePileup(f, snpChan)
	}()

	piArr := []Pi{}
	for snp := range snpChan {
		if filterSNP(snp) {
			pi := CalcPi(snp)
			if !math.IsNaN(pi) {
				piArr = append(piArr, pi)
			}
		}
	}
	fmt.Println("Finish calculating pi!")

	totalPis := divivePis(piArr, 100)
	means, vars, ns := collectCov(totalPis, profile, parsePos(pos), maxl)
	ks, ksvar := collectKs(totalPis, profile, parsePos(pos))

	w, err := os.Create(outFile)
	if err != nil {
		panic(err)
	}
	defer w.Close()

	w.WriteString(fmt.Sprintf("# ks = %g, ksvar = %g\n", ks, ksvar))

	for i := 0; i < len(means); i++ {
		w.WriteString(fmt.Sprintf("%d\t%g\t%g\t%d\n", i, means[i], vars[i], ns[i]))
	}
}

func filterSNP(snp SNP) (good bool) {
	num := 0
	for i := 0; i < len(snp.ReadBases); i++ {
		if snp.ReadBases[i] != '*' {
			num++
		}
	}
	if num >= 10 {
		good = true
	}

	return
}

func parsePos(pos int) byte {
	var p byte
	switch pos {
	case 0:
		p = NonCoding
	case 1:
		p = FirstPos
	case 2:
		p = SecondPos
	case 3:
		p = ThirdPos
	case 4:
		p = FourFold
	default:
		p = ThirdPos
	}

	return p
}

func generatePis(snaArr []SNP) []Pi {
	piArr := []Pi{}
	for _, snp := range snaArr {
		piArr = append(piArr, CalcPi(snp))
	}
	return piArr
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

func collectCov(totalPis [][]Pi, profile []byte, pos byte, maxl int) (means, vars []float64, ns []int) {
	mvs := make([]*meanvar.MeanVar, maxl)
	for i := 0; i < len(mvs); i++ {
		mvs[i] = meanvar.New()
	}

	for i := 0; i < len(totalPis); i++ {
		corrs, ns := CalcCovRate(totalPis[i], profile, pos, maxl)
		for l := 0; l < maxl; l++ {
			if ns[l] >= 10 && !math.IsNaN(corrs[l]) {
				mvs[l].Increment(corrs[l])
			}
		}
	}

	means = make([]float64, len(mvs))
	vars = make([]float64, len(mvs))
	ns = make([]int, len(mvs))
	for i := 0; i < len(mvs); i++ {
		if !math.IsNaN(mvs[i].Var.GetResult()) {
			means[i] = mvs[i].Mean.GetResult()
			ns[i] = mvs[i].Mean.GetN()
			vars[i] = mvs[i].Var.GetResult()
		}
	}

	return
}

func collectKs(totalPis [][]Pi, profile []byte, pos byte) (m, v float64) {
	mv := meanvar.New()
	for _, piArr := range totalPis {
		ks := CalcKs(piArr, profile, pos)
		mv.Increment(ks)
	}
	m = mv.Mean.GetResult()
	v = mv.Var.GetResult()
	return
}
