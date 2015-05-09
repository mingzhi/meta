package cov

import (
	"github.com/mingzhi/gomath/stat/correlation"
	"github.com/mingzhi/gomath/stat/desc"
)

type CovCalculator struct {
	corrs []*correlation.BivariateCovariance
}

func NewCovCalculator(maxl int, bias bool) *CovCalculator {
	cc := CovCalculator{}
	cc.corrs = make([]*correlation.BivariateCovariance, maxl)
	for i := 0; i < maxl; i++ {
		cc.corrs[i] = correlation.NewBivariateCovariance(bias)
	}
	return &cc
}

func (cc *CovCalculator) Increment(i int, x, y float64) {
	cc.corrs[i].Increment(x, y)
}

func (cc *CovCalculator) GetResult(i int) float64 {
	return cc.corrs[i].GetResult()
}

func (cc *CovCalculator) GetMeanXY(i int) float64 {
	return cc.corrs[i].MeanX() * cc.corrs[i].MeanY()
}

func (cc *CovCalculator) GetN(i int) int {
	return cc.corrs[i].GetN()
}

func (cc *CovCalculator) Append(cc2 *CovCalculator) {
	for i := 0; i < len(cc.corrs); i++ {
		cc.corrs[i].Append(cc2.corrs[i])
	}
}

type KsCalculator struct {
	Mean *desc.Mean
	Var  *desc.Variance
}

func NewKsCalculator() *KsCalculator {
	kc := KsCalculator{}
	kc.Mean = desc.NewMean()
	kc.Var = desc.NewVarianceWithBiasCorrection()
	return &kc
}

func (kc *KsCalculator) Increment(v float64) {
	kc.Mean.Increment(v)
	kc.Var.Increment(v)
}

func (kc *KsCalculator) Append(kc2 *KsCalculator) {
	kc.Mean.Append(kc2.Mean)
	kc.Var.Append(kc2.Var)
}
