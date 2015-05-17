package cov

import (
	"github.com/mingzhi/meta/genome"
	"math"
)

const (
	AlphabetDNA = "ATGCatgc"
)

// Generate substitution profile according to the position profile.
func SubProfile(read, nucl, profile []byte, pos int) []float64 {
	// determine the codon position.
	var cp byte
	switch pos {
	case 1:
		cp = genome.FirstPos
	case 2:
		cp = genome.SecondPos
	case 3:
		cp = genome.ThirdPos
	case 4:
		cp = genome.FourFold
	}
	subs := make([]float64, len(nucl))
	for i := 0; i < len(subs); i++ {
		p := profile[i]
		var match bool
		if cp == genome.ThirdPos {
			match = p == genome.ThirdPos || p == genome.FourFold
		} else {
			match = p == cp
		}
		valid := validate(read[i]) && validate(nucl[i])
		if match && valid {
			if read[i] == nucl[i] {
				subs[i] = 0
			} else {
				subs[i] = 1
			}
		} else {
			subs[i] = math.NaN()
		}
	}
	return subs
}

func validate(r byte) bool {
	for i := 0; i < len(AlphabetDNA); i++ {
		if AlphabetDNA[i] == r {
			return true
		}
	}
	return false
}
