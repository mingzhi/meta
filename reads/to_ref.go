package reads

import (
	"bytes"
	"github.com/biogo/hts/sam"
)

// Obtain the sequence of a read mapping to the reference genome.
// Return the mapped sequence.
func Map2Ref(r *sam.Record) []byte {
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
func MapMated2Ref(r PairedEndRead) []byte {
	s1 := Map2Ref(r.ReadLeft)
	s2 := Map2Ref(r.ReadRight)
	space := r.ReadRight.Pos - (r.ReadLeft.Pos + len(s1))
	if space > 0 {
		s1 = append(s1, bytes.Repeat([]byte{'*'}, space)...)
		s1 = append(s1, s2...)
	} else {
		s1 = append(s1, s2[-space:]...)
	}

	return s1
}
