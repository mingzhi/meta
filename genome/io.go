package genome

import (
	"github.com/mingzhi/biogo/seq"
	"io/ioutil"
	"os"
	"path/filepath"
)

func readProfile(fileName string) Profile {
	f, err := os.Open(fileName)
	if err != nil {
		panic(err)
	}
	defer f.Close()

	data, err := ioutil.ReadAll(f)
	if err != nil {
		panic(err)
	}

	return Profile(data)
}

func readFasta(fileName string) []byte {
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

	return seqs[0].Seq
}

// Load position profile to the genome,
// from the file in base folder.
func LoadProfile(g *Genome, base string) {
	fileName := filepath.Join(base, g.RefAcc()+".pos")
	g.PosProfile = readProfile(fileName)
}

// Load sequence to the genome,
// from the .fna file in base folder.
func LoadFna(g *Genome, base string) {
	fileName := filepath.Join(base, g.RefAcc()+".fna")
	g.Seq = readFasta(fileName)
}
