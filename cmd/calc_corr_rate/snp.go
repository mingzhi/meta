package main

// SNP contains single nucleotide polymorphisms.
type SNP struct {
	Genome    string
	Position  int // 1-coordinate system [1 - N]
	RefBase   byte
	Number    int
	ReadBases []byte
	BaseQuals []int
}
