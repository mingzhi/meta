package cov

import (
	"github.com/mingzhi/gomath/stat/desc"
	"math"
)

func getPosIndices(subs []float64) (ints []int) {
	for i := 0; i < len(subs); i++ {
		if !math.IsNaN(subs[i]) {
			ints = append(ints, i)
		}
	}

	return
}

// Calculate correlations of subsitutions.
func SubCorr(subs []float64, cc *CovCalculator, kc *KsCalculator, maxl int) {
	// Calculate kc and obtain position indices.
	ints := getPosIndices(subs)
	for _, i := range ints {
		kc.Increment(subs[i])
		cc.Increment(0, subs[i], subs[i])
	}

	// Calculate cc.
	for j := 0; j < len(ints); j++ {
		for k := j + 1; k < len(ints); k++ {
			l := ints[k] - ints[j]
			if l >= maxl {
				break
			} else {
				x, y := subs[ints[j]], subs[ints[k]]
				cc.Increment(l, x, y)
			}
		}
	}
}

// Calculate structure covariance for a subsitution matrix.
func SubMatrixSCov(subMatrix [][]float64, cs *MeanCovCalculator, maxl int) {
	ints := getPosIndices(subMatrix[0])
	for i := 0; i < len(ints); i++ {
		for j := i; j < len(ints); j++ {
			l := ints[j] - ints[i]
			xs, ys := []float64{}, []float64{}
			for k := 0; k < len(subMatrix); k++ {
				x, y := subMatrix[k][i], subMatrix[k][j]
				xs = append(xs, x)
				ys = append(ys, y)
			}
			if len(xs) > 3 {
				cs.Increment(xs, ys, l)
			}
		}
	}
}

// Calculate mutation covariance for a substitution matrix.
func SubMatrixMCov(subMatrix [][]float64, cm *MeanCovCalculator, maxl int) {
	ints := getPosIndices(subMatrix[0])

	for k := 0; k < len(subMatrix); k++ {
		m := make(map[int][][]float64)
		for i := 0; i < len(ints); i++ {
			for j := i; j < len(ints); j++ {
				l := ints[j] - ints[i]
				if _, found := m[l]; !found {
					m[l] = make([][]float64, 2)
				}
				x, y := subMatrix[k][i], subMatrix[k][j]
				m[l][0] = append(m[l][0], x)
				m[l][1] = append(m[l][1], y)
			}
		}

		for l, values := range m {
			if len(values[0]) > 3 {
				cm.Increment(values[0], values[1], l)
			}
		}
	}
}

// Calculate rate mutation covariance for a substituion matrix.
func SubMatrixRCov(subMatrix [][]float64, cr *MeanCovCalculator, maxl int) {
	means := make([]*desc.Mean, len(subMatrix[0]))
	for i := 0; i < len(means); i++ {
		means[i] = desc.NewMean()
	}

	for i := 0; i < len(subMatrix); i++ {
		for j := 0; j < len(subMatrix[i]); j++ {
			means[j].Increment(subMatrix[i][j])
		}
	}

	subs := []float64{}
	for i := 0; i < len(means); i++ {
		subs = append(subs, means[i].GetResult())
	}

	SubMatrixMCov([][]float64{subs}, cr, maxl)
}
