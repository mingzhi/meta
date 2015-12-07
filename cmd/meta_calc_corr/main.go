package main

import (
	"flag"
	"fmt"
	"github.com/mingzhi/ncbiftp/taxonomy"
	"os"
)

var (
	bamFileName  string
	genomeAcc    string
	outFile      string
	maxl         int
	codonTableId string
)

func init() {
	flag.StringVar(&codonTableId, "codon", "11", "codon table id")
	flag.IntVar(&maxl, "maxl", 500, "maxl")
	flag.Parse()
	bamFileName = flag.Arg(0)
	genomeAcc = flag.Arg(1)
	outFile = flag.Arg(2)
}

func main() {
	codonTable := taxonomy.GeneticCodes()[codonTableId]
	genomeFile := genomeAcc + ".fna"
	pttFile := genomeAcc + ".ptt"
	profile := ProfileGenome(genomeFile, pttFile, codonTable)
	fmt.Println("Finish generating profile!")
	_, samRecordChan := ReadBamFile(bamFileName)
	snpChan := Pileup(samRecordChan)
	c := Calc(snpChan, profile, FourFold, maxl)
	means, covs := Collect(c)

	w, err := os.Create(outFile)
	if err != nil {
		panic(err)
	}
	defer w.Close()

	for i := 0; i < len(means); i++ {
		w.WriteString(fmt.Sprintf("%d\t%g\t%g\n", i, means[i], covs[i]))
	}
}
