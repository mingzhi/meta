package cov

import (
	"github.com/mingzhi/gomath/stat/desc"
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

// Calculate structure covariance for a subsitution matrix.
func SubMatrixSCov(subMatrix [][]float64, cs *MeanCovCalculator, maxl int) {
	ncol := len(subMatrix[0])
	nrow := len(subMatrix)
	for l := 0; l < maxl; l++ {
		for i := 0; i < ncol-l; i++ {
			xs, ys := []float64{}, []float64{}
			for k := 0; k < nrow; k++ {
				x, y := subMatrix[k][i], subMatrix[k][i+l]
				if !math.IsNaN(x) && !math.IsNaN(y) {
					xs = append(xs, x)
					ys = append(ys, y)
				}
			}
			if len(xs) > 3 {
				cs.Increment(xs, ys, l)
			}
		}
	}
}

// Calculate mutation covariance for a substitution matrix.
func SubMatrixMCov(subMatrix [][]float64, cm *MeanCovCalculator, maxl int) {
	for i := 0; i < len(subMatrix); i++ {
		for l := 0; l < maxl; l++ {
			xs, ys := []float64{}, []float64{}
			for k := 0; k < len(subMatrix[i])-l; k++ {
				x, y := subMatrix[i][k], subMatrix[i][k+l]
				if !math.IsNaN(x) && !math.IsNaN(y) {
					xs = append(xs, x)
					ys = append(ys, y)
				}
			}
			if len(xs) > 3 {
				cm.Increment(xs, ys, l)
			}
		}
	}
}

// Calculate rate mutation covariance for a substituion matrix.
func SubMatrixRCov(subMatrix [][]float64, cr *MeanCovCalculator, maxl int) {
	means := make([]*desc.Mean, len(subMatrix[0]))
	for i := 0; i < len(subMatrix); i++ {
		for j := 0; j < len(subMatrix[i]); j++ {
			means[j].Increment(subMatrix[i][j])
		}
	}

	for l := 0; l < maxl; l++ {
		xs, ys := []float64{}, []float64{}
		for k := 0; k < len(means)-l; k++ {
			x, y := means[k].GetResult(), means[k+l].GetResult()
			if !math.IsNaN(x) && !math.IsNaN(y) {
				xs = append(xs, x)
				ys = append(ys, y)
			}
		}
		if len(xs) > 3 {
			cr.Increment(xs, ys, l)
		}
	}
}
