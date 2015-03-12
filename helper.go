package meta

import (
	"bytes"
	"code.google.com/p/biogo.hts/sam"
	"github.com/mingzhi/biogo/seq"

	"io/ioutil"
	"os"
)

// Helper functions

// Obtain the sequence of a read mapping to the reference genome.
// Return the mapped sequence.
func map2Ref(r *sam.Record) []byte {
	s := []byte{}
	p := 0                 // position in the read sequence.
	read := r.Seq.Expand() // read sequence.
	for _, c := range r.Cigar {
		switch c.Type() {
		case sam.CigarMatch, sam.CigarMismatch, sam.CigarEqual:
			s = append(s, read[p:p+c.Len()]...)
			p += c.Len()
		case sam.CigarInsertion, sam.CigarSoftClipped, sam.CigarHardClipped:
			p += c.Len()
		case sam.CigarDeletion, sam.CigarSkipped:
			s = append(s, bytes.Repeat([]byte{'*'}, c.Len())...)
		}
	}

	return s
}

// Obtain the sequence of a pair of mated reads to the reference genome.
func mated2Ref(r PairedEndRead) []byte {
	s1 := map2Ref(r.ReadLeft)
	s2 := map2Ref(r.ReadRight)
	space := r.ReadRight.Pos - (r.ReadLeft.Pos + len(s1))
	if space > 0 {
		s1 = append(s1, bytes.Repeat([]byte{'*'}, space)...)
		s1 = append(s1, s2...)
	} else {
		s1 = append(s1, s2[-space:]...)
	}

	return s1
}

// Read postion profile for a genome.
func ReadPosProfile(fileName string) []byte {
	f, err := os.Open(fileName)
	if err != nil {
		panic(err)
	}
	defer f.Close()
	data, err := ioutil.ReadAll(f)
	if err != nil {
		panic(err)
	}

	return data
}

// Read fasta sequence
func ReadFasta(fileName string) *seq.Sequence {
	f, err := os.Open(fileName)
	if err != nil {
		panic(err)
	}
	defer f.Close()

	rd := seq.NewFastaReader(f)

	seqs, err := rd.ReadAll()
	if err != nil {
		panic(err)
	}

	return seqs[0]
}

// return max int
func maxInt(a, b int) int {
	if a > b {
		return a
	} else {
		return b
	}
}

// return min int
func minInt(a, b int) int {
	if a < b {
		return a
	} else {
		return b
	}
}
