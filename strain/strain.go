package strain

import (
	"github.com/mingzhi/biogo/seq"
	"github.com/mingzhi/meta/genome"
	"github.com/mingzhi/ncbiftp/seqrecord"
	"github.com/mingzhi/ncbiftp/taxonomy"
	"log"
	"os"
	"path/filepath"
	"strings"
)

// Strain information.
type Strain struct {
	Name        string          // taxonomy name.
	TaxId       string          // taxonomy Id.
	ProjectId   string          // project Id.
	Genomes     []genome.Genome // genomes.
	Path        string          // path in referenc genome diretory.
	GeneticCode string          // genetic code id.
	Species     string          // species name
	Status      string          // status of sequencing.
}

// Position profile genomes.
// base is where the genome folder in NCBI ftp.
func (s *Strain) ProfileGenomes(base string, gcMap map[string]*taxonomy.GeneticCode) {
	// absolute path storing the genome.
	dir := filepath.Join(base, s.Path)
	// genetic codon table to determine four-fold sites.
	gc := gcMap[s.GeneticCode]
	// for each genome.
	for _, g := range s.Genomes {
		// read genome sequence.
		fnaFileName := g.RefAcc() + ".fna"
		fnaFilePath := filepath.Join(dir, fnaFileName)
		genomeSeq := readFasta(fnaFilePath)[0].Seq

		// initialize position profile.
		profile := make(genome.Profile, len(genomeSeq))

		// read protein features from ptt file.
		pttFileName := g.RefAcc() + ".ptt"
		pttFilePath := filepath.Join(dir, pttFileName)
		pttFile := seqrecord.NewPttFile(pttFilePath)
		ptts := pttFile.ReadAll()

		// for each gene (a ptt), record codon position.
		for _, ptt := range ptts {
			// Prepare nucleotide sequence,
			// we need it for determine 4-fold codons.
			var nucl []byte
			if ptt.Loc.To >= ptt.Loc.From {
				nucl = genomeSeq[ptt.Loc.From-1 : ptt.Loc.To]
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
					prof[j] = genome.FirstPos
				case 2:
					prof[j] = genome.SecondPos
				case 0:
					codon := nucl[j-2 : j+1]
					if gc.FFCodons[string(codon)] {
						prof[j] = genome.FourFold
					} else {
						prof[j] = genome.ThirdPos
					}
				}
			}

			if ptt.Loc.Strand == "-" {
				prof = seq.Reverse(prof)
			}

			for j, p := range prof {
				profile[ptt.Loc.From-1+j] = p
			}
		}

		// Save genome profile to file.
		fileName := strings.Replace(fnaFilePath, "fna", "pos", -1)
		writePosProfile(fileName, profile)
	}
}

func readFasta(fileName string) []*seq.Sequence {
	f, err := os.Open(fileName)
	if err != nil {
		log.Panic(err)
	}
	defer f.Close()

	fastaRd := seq.NewFastaReader(f)
	seqs, err := fastaRd.ReadAll()
	if err != nil {
		log.Panic(err)
	}

	return seqs
}

func writePosProfile(fileName string, profile []byte) {
	w, err := os.Create(fileName)
	if err != nil {
		log.Panic(err)
	}
	defer w.Close()

	w.Write(profile)
}
