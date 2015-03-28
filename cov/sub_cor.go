package cov

import (
	"math"
)

// Calculate correlations of subsitutions.
func SubCorr(subs []float64, cc *CovCalculator, kc *KsCalculator, maxl int) {
	// Calculate kc and obtain position indices.
	posIndices := []int{}
	for k := 0; k < len(subs); k++ {
		if !math.IsNaN(subs[k]) {
			kc.Increment(subs[k])
			posIndices = append(posIndices, k)
			cc.Increment(0, subs[k], subs[k])
		}
	}
	// Calculate cc.
	for j := 0; j < len(posIndices); j++ {
		for k := j + 1; k < len(posIndices); k++ {
			l := posIndices[k] - posIndices[j]
			if l >= maxl {
				break
			} else {
				x, y := subs[posIndices[j]], subs[posIndices[k]]
				cc.Increment(l, x, y)
			}
		}
	}
}
