package main

import (
	"github.com/mingzhi/gomath/stat/correlation"
)

// Calculate covariance of rates.
func CalcCovRate(pis []Pi, profile []byte, pos byte, maxl int) (covs []float64, n []int) {
	corrs := make([]*correlation.BivariateCovariance, maxl)
	for i := 0; i < maxl; i++ {
		corrs[i] = correlation.NewBivariateCovariance(false)
	}

	for i := 0; i < len(pis); i++ {
		pos1 := profile[pis[i].Position-1]
		pos1 = convertPos(pos1, pos)

		if pos1 == pos {
			for j := i; j < len(pis); j++ {
				pos2 := profile[pis[j].Position-1]
				pos2 = convertPos(pos2, pos)

				if pos2 == pos {
					distance := pis[j].Position - pis[i].Position
					if distance < maxl {
						corrs[distance].Increment(pis[i].Pi, pis[j].Pi)
					} else {
						break
					}
				}
			}
		}

	}

	for i := 0; i < maxl; i++ {
		covs = append(covs, corrs[i].GetResult())
		n = append(n, corrs[i].GetN())
	}

	return
}

func convertPos(pos1, pos byte) byte {
	if pos == ThirdPos && pos1 == FourFold {
		pos1 = ThirdPos
	}

	return pos1
}
