package main

type SNP struct {
	Genome    string
	Position  int
	RefBase   byte
	Number    int
	ReadBases []byte
	BaseQuals []int
}
