// This program calculate rate correlations (c_R) and structure correlation (c_R),
// from a mapping results of metagenomic sequences to a reference genome.
// We need two inputs:
// 1. the mapping results in .bam format, sorted;
// 2. the reference genome sequence
// 3. the protein features.
package main

import (
	"flag"
	"fmt"
	"github.com/mingzhi/ncbiftp/genomes/profiling"
	"github.com/mingzhi/ncbiftp/taxonomy"
	"log"
	"math"
	"os"
	"runtime/pprof"
)

var (
	bamFileName  string // mapping results in .bam file
	genomeFile   string // genome accession number
	pttFile      string // protein feature file
	outFile      string // output file.
	maxl         int    // max length of correlation.
	codonTableID string // codon table ID.
	pos          int    // position for calculation.
	minBQ        int
	minDepth     int
	minMQ        int
	cpuprofile   string
)

func init() {
	flag.StringVar(&codonTableID, "codon", "11", "Codon table ID")
	flag.IntVar(&maxl, "maxl", 500, "Maximum length of correlation distance")
	flag.IntVar(&pos, "pos", 4, "Position for SNP calculation")
	flag.IntVar(&minBQ, "min-BQ", 13, "Minimum base quality for a base to be considered")
	flag.IntVar(&minDepth, "min-depth", 20, "At a position, mimimum number of reads included to calculation")
	flag.IntVar(&minMQ, "min-MQ", 0, "Minimum read mapping quality")
	flag.StringVar(&cpuprofile, "cpuprofile", "", "write cpu profile to file")
	flag.Parse()
	if flag.NArg() < 4 {
		fmt.Println("meta_calc_corr <bam file> <ref genome sequence> <protein feature file> <output file>")
		os.Exit(1)
	}
	bamFileName = flag.Arg(0)
	genomeFile = flag.Arg(1)
	pttFile = flag.Arg(2)
	outFile = flag.Arg(3)
}

func main() {
	// Profile CPU usage.
	if cpuprofile != "" {
		f, err := os.Create(cpuprofile)
		if err != nil {
			log.Fatal(err)
		}
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	}

	// Obtain codon table for identifying four-fold degenerate sites.
	codonTable := taxonomy.GeneticCodes()[codonTableID]
	// Profiling genome using reference sequence and protein feature data.
	profile := profiling.ProfileGenome(genomeFile, pttFile, codonTable)

	// Read mapping records in sam formate from the .bam file.
	_, samRecordChan := ReadBamFile(bamFileName)
	// Pileup the mapped bases for each genomic position.
	snpChan := Pileup(samRecordChan)
	// Using the pileup data for correlation calculation.
	positionType := convertPosType(pos)
	cChan := Calc(snpChan, profile, positionType, maxl)
	// Collect results from the calculator.
	cSs, cRs, xBars, yBars, cTs := Collect(maxl, cChan)

	w, err := os.Create(outFile)
	if err != nil {
		panic(err)
	}
	defer w.Close()

	// Write header.
	w.WriteString("#i\tcs\tvar cs\tn cs\tcr\tvar cr\tn cr\tct\tvar ct\tn ct\txbar\tvar xbar\tn xbar\tybar\tvar ybar\tn ybar\n")
	for i := 0; i < len(cSs); i++ {
		f := cRs[i].Var.GetResult()
		n := cRs[i].Var.GetN()
		if !math.IsNaN(f) && n > 0 {
			w.WriteString(
				fmt.Sprintf("%d\t%g\t%g\t%d\t%g\t%g\t%d\t%g\t%g\t%d\t%g\t%g\t%d\t%g\t%g\t%d\n",
					i,
					cSs[i].Mean.GetResult(), cSs[i].Var.GetResult(), cSs[i].Mean.GetN(),
					cRs[i].Mean.GetResult(), cRs[i].Var.GetResult(), cRs[i].Mean.GetN(),
					cTs[i].Mean.GetResult(), cTs[i].Var.GetResult(), cTs[i].Mean.GetN(),
					xBars[i].Mean.GetResult(), xBars[i].Var.GetResult(), xBars[i].Mean.GetN(),
					yBars[i].Mean.GetResult(), yBars[i].Var.GetResult(), yBars[i].Mean.GetN(),
				),
			)
		}
	}
}

func convertPosType(pos int) byte {
	var p byte
	switch pos {
	case 0:
		p = profiling.NonCoding
	case 1:
		p = profiling.FirstPos
		break
	case 2:
		p = profiling.SecondPos
		break
	case 3:
		p = profiling.ThirdPos
		break
	case 4:
		p = profiling.FourFold
		break
	default:
		p = profiling.Coding
	}

	return p
}
