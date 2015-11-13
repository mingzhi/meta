package main

import (
	"bytes"
	"github.com/mingzhi/gomath/stat/desc"
)

type Pi struct {
	Position int
	Pi       float64
}

// Calc Pi from each snp.
func CalcPi(snp SNP) (pi Pi) {
	pi.Position = snp.Position
	pi.Pi = calcPi(snp.ReadBases)
	return pi
}

func calcPi(bases []byte) (pi float64) {
	// convert bases to upper case.
	upperBases := bytes.ToUpper(bases)

	mean := desc.NewMean()
	for i := 0; i < len(upperBases); i++ {
		if upperBases[i] == '*' {
			continue
		}

		for j := i + 1; j < len(upperBases); j++ {
			if upperBases[j] == '*' {
				continue
			}

			if upperBases[i] != upperBases[j] {
				mean.Increment(1.0)
			} else {
				mean.Increment(0.0)
			}
		}
	}

	pi = mean.GetResult()

	return
}
