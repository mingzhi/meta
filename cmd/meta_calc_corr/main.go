// This program calculate rate correlations (c_R) and structure correlation (c_R),
// from a mapping results of metagenomic sequences to a reference genome.
// We need two inputs:
// 1. the mapping results in .bam format;
// 2. the reference genome sequence and the protein features.
package main

import (
	"flag"
	"fmt"
	"github.com/mingzhi/ncbiftp/taxonomy"
	"os"
)

var (
	bamFileName        string // mapping results in .bam file
	genomeFile         string // genome accession number
	proteinFeatureFile string // protein feature file
	outFile            string
	maxl               int
	codonTableId       string
	pos                int // position for calculation.
)

func init() {
	flag.StringVar(&codonTableId, "codon", "11", "codon table id")
	flag.IntVar(&maxl, "maxl", 500, "maxl")
	flag.IntVar(&pos, "pos", 4, "position for SNP calculation.")
	flag.Parse()
	if flag.NArg() < 4 {
		fmt.Println("meta_calc_corr <bam file> <ref genome sequence> <protein feature file> <output file>")
		os.Exit(1)
	}
	bamFileName = flag.Arg(0)
	genomeFile = flag.Arg(1)
	proteinFeatureFile = flag.Arg(2)
	outFile = flag.Arg(3)
}

func main() {
	// Obtain codon table for following genome profiling.
	codonTable := taxonomy.GeneticCodes()[codonTableId]
	// Profiling genome using reference sequence and protein feature data.
	profile := ProfileGenome(genomeFile, proteinFeatureFile, codonTable)

	// Read mapping records in sam formate from the .bam file.
	_, samRecordChan := ReadBamFile(bamFileName)
	// Pileup the mapped bases for each genomic position.
	snpChan := Pileup(samRecordChan)
	// Using the pileup data for correlation calculation.
	geneLength := len(profile) / 1000
	cChan := Calc(snpChan, profile, convertPos(pos), maxl, geneLength)
	// Collect results from the calculator.
	means, covs := Collect(maxl, cChan)

	w, err := os.Create(outFile)
	if err != nil {
		panic(err)
	}
	defer w.Close()

	for i := 0; i < len(means); i++ {
		w.WriteString(fmt.Sprintf("%d\t%g\t%g\t%d\t%g\t%g\t%d\n", i, means[i].Mean.GetResult(), means[i].Var.GetResult(), means[i].Mean.GetN(),
			covs[i].Mean.GetResult(), covs[i].Var.GetResult(), covs[i].Mean.GetN()))
	}
}

func convertPos(pos int) byte {
	var p byte
	switch pos {
	case 1:
		p = FirstPos
	case 2:
		p = SecondPos
	case 3:
		p = ThirdPos
	case 4:
		p = FourFold
	}

	return p
}
