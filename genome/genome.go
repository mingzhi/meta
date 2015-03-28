package genome

import (
	"regexp"
)

type Genome struct {
	Accession  string  // RefSeq accession.
	Replicon   string  // type of DNA replication.
	Length     int     // length of the genome.
	Seq        []byte  // genome sequence.
	PosProfile Profile // position profile.
}

type Profile []byte

const (
	FirstPos byte = 1 << iota
	SecondPos
	ThirdPos
	FourFold
)

func (g Genome) RefAcc() string {
	return FindRefAcc(g.Accession)
}

func FindRefAcc(s string) string {
	r := regexp.MustCompile("\\w\\w_\\w+\\d+")
	return r.FindString(s)
}
