package main

import (
	"github.com/mingzhi/gomath/stat/desc"
)

func CalcKs(pis []Pi, profile []byte, pos byte) float64 {
	mean := desc.NewMean()
	for i := 0; i < len(pis); i++ {
		pos1 := profile[pis[i].Position]
		pos1 = convertPos(pos1, pos)
		if pos1 == pos {
			mean.Increment(pis[i].Pi)
		}

	}
	return mean.GetResult()
}
