package main

import (
	"sort"

	"github.com/mingzhi/ncbiftp/taxonomy"
)

// Codon stores a codon value, the position in a genome and the read id.
type Codon struct {
	Seq     string
	ReadID  int
	GenePos int
}

// CodonPile stores a pile of Codon, which are at a particular genome position.
type CodonPile struct {
	Codons []Codon
}

// NewCodonPile return a new CodonPile.
func NewCodonPile() *CodonPile {
	return &CodonPile{}
}

// Append appends a new Codon.
func (cp *CodonPile) Append(c Codon) {
	cp.Codons = append(cp.Codons, c)
}

// SortReadByID sorts Codons by ReadName
func (cp *CodonPile) SortReadByID() {
	sort.Slice(cp.Codons, func(i, j int) bool { return cp.Codons[i].ReadID < cp.Codons[j].ReadID })
}

// SearchReadByID search a codon by ReadName. If not found, it returns nil.
func (cp *CodonPile) SearchReadByID(readID int) Codon {
	data := cp.Codons
	i := sort.Search(len(data), func(i int) bool { return data[i].ReadID >= readID })
	if i < len(data) && data[i].ReadID == readID {
		return data[i]
	}
	return Codon{ReadID: -1}
}

// Len return the lenght of pileup Codons.
func (cp *CodonPile) Len() int {
	return len(cp.Codons)
}

// CodonGene represents a gene with an array of CodonPile.
type CodonGene struct {
	CodonPiles []*CodonPile
}

// NewCodonGene return a new CodonGene.
func NewCodonGene() *CodonGene {
	return &CodonGene{}
}

// AddCodon add a codon.
func (cg *CodonGene) AddCodon(c Codon) {
	for len(cg.CodonPiles) <= c.GenePos {
		cg.CodonPiles = append(cg.CodonPiles, NewCodonPile())
	}
	cg.CodonPiles[c.GenePos].Append(c)
}

// DepthAt return the pile depth at position i.
func (cg *CodonGene) DepthAt(i int) int {
	if len(cg.CodonPiles) <= i {
		return 0
	}
	return cg.CodonPiles[i].Len()
}

// Len returns length of CodonPile array.
func (cg *CodonGene) Len() int {
	return len(cg.CodonPiles)
}

// SortCodonByReadID sorts codons by read id for each codon pile.
func (cg *CodonGene) SortCodonByReadID() {
	for i := range cg.CodonPiles {
		cg.CodonPiles[i].SortReadByID()
	}
}

// CodonPair stores a pair of Codon
type CodonPair struct {
	A, B Codon
}

// PairCodonAt pairs codons at positions i and j.
func (cg *CodonGene) PairCodonAt(i, j int) (pairs []CodonPair) {
	if i >= len(cg.CodonPiles) || j >= len(cg.CodonPiles) {
		return
	}

	pile1 := cg.CodonPiles[i]
	pile2 := cg.CodonPiles[j]
	for k := 0; k < pile1.Len(); k++ {
		codon1 := pile1.Codons[k]
		codon2 := pile2.SearchReadByID(codon1.ReadID)

		if codon2.ReadID != -1 {
			pairs = append(pairs, CodonPair{A: codon1, B: codon2})
		}
		if len(pairs) == pile2.Len() {
			break
		}
	}
	return
}

// SynoumousSplitCodonPairs split codon pairs into synoumous pairs.
func SynoumousSplitCodonPairs(codonPairs []CodonPair, codeTable *taxonomy.GeneticCode) [][]CodonPair {
	var splittedPairs [][]CodonPair
	var aaArray []string
	for _, codonPair := range codonPairs {
		hasGap := false
		for _, codon := range []Codon{codonPair.A, codonPair.B} {
			for _, b := range codon.Seq {
				if !isATGC(byte(b)) {
					hasGap = true
					break
				}
			}
			if hasGap {
				break
			}
		}

		if hasGap {
			continue
		}

		a := codeTable.Table[codonPair.A.Seq]
		b := codeTable.Table[codonPair.B.Seq]
		ab := string([]byte{a, b})
		index := -1
		for i, aa := range aaArray {
			if aa == ab {
				index = i
			}
		}
		if index == -1 {
			index = len(aaArray)
			aaArray = append(aaArray, ab)
			splittedPairs = append(splittedPairs, []CodonPair{})
		}
		splittedPairs[index] = append(splittedPairs[index], codonPair)
	}
	return splittedPairs
}
