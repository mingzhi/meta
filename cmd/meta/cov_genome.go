package main

import (
	"flag"
	"fmt"
	"github.com/mingzhi/gomath/stat/desc"
	"github.com/mingzhi/meta"
	"log"
	"math"
	"os"
	"path/filepath"
	"strings"
)

type cmdCovGenome struct {
	workspace *string // workspace
	prefix    *string // prefix of bacterial species
	ref       *string // reference genome database.
	samout    *string // sam output diretory
	maxl      *int    // max length of correlation
}

func (cmd *cmdCovGenome) Flags(fs *flag.FlagSet) *flag.FlagSet {
	cmd.workspace = fs.String("w", "", "workspace")
	cmd.ref = fs.String("r", "", "reference genome diretory")
	cmd.samout = fs.String("s", "", "sam output folder")
	cmd.prefix = fs.String("p", "", "prefix")
	cmd.maxl = fs.Int("m", 100, "maxl")
	return fs
}

func (cmd *cmdCovGenome) Run(args []string) {
	registerLogger()
	speciesMapFileName := filepath.Join(*cmd.workspace, "species_map.json")
	speciesMap := meta.ReadSpeciesMap(speciesMapFileName)
	strains, found := speciesMap[*cmd.prefix]
	if !found {
		log.Fatalf("Can not find strain information for %s from species_map.json\n", *cmd.prefix)
	}

	var cm []*MeanVar
	for i := 0; i < *cmd.maxl; i++ {
		cm = append(cm, NewMeanVar())
	}

	for _, s := range strains {
		for _, genome := range s.Genomes {
			if strings.Contains(genome.Replicon, "chromosome") {
				bamFileName := s.Path + ".align.bam"
				bamFilePath := filepath.Join(*cmd.samout, s.Path, bamFileName)
				_, records := meta.ReadBamFile(bamFilePath)
				// Read postion profile for the genome.
				posFileName := findRefAcc(genome.Accession) + ".pos"
				posFilePath := filepath.Join(*cmd.ref, s.Path, posFileName)
				genome.PosProfile = meta.ReadPosProfile(posFilePath)
				// Read sequence for the genome.
				fnaFileName := findRefAcc(genome.Accession) + ".fna"
				fnaFilePath := filepath.Join(*cmd.ref, s.Path, fnaFileName)
				genome.Seq = meta.ReadFasta(fnaFilePath).Seq

				cc := meta.CovGenome(records, genome, *cmd.maxl)
				for i := 0; i < *cmd.maxl; i++ {
					if !math.IsNaN(cc.GetResult(i)) {
						cm[i].Increment(cc.GetResult(i))
					}

				}

				outFileName := s.Path + "_cov.tsv"
				outFilePath := filepath.Join(*cmd.workspace, outFileName)
				out, err := os.Create(outFilePath)
				if err != nil {
					panic(err)
				}
				defer out.Close()

				for i := 0; 3*i < *cmd.maxl; i++ {
					out.WriteString(fmt.Sprintf("%d\t%f\t%d\n", i, cc.GetResult(3*i), cc.GetN(3*i)))
				}
			}

		}
	}

	outFileName := *cmd.prefix + "_cov.tsv"
	outFilePath := filepath.Join(*cmd.workspace, outFileName)
	out, err := os.Create(outFilePath)
	if err != nil {
		panic(err)
	}
	defer out.Close()

	for i := 0; 3*i < *cmd.maxl; i++ {
		out.WriteString(fmt.Sprintf("%d\t%f\t%d\n", i, cm[3*i].Mean.GetResult(), cm[3*i].Mean.GetN()))
	}
}

type MeanVar struct {
	Mean *desc.Mean
	Var  *desc.Variance
}

func NewMeanVar() *MeanVar {
	mean := desc.NewMean()
	variance := desc.NewVarianceWithBiasCorrection()
	return &MeanVar{
		Mean: mean,
		Var:  variance,
	}
}

func (m *MeanVar) Increment(d float64) {
	m.Mean.Increment(d)
	m.Var.Increment(d)
}
