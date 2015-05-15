package cov

import (
	"github.com/mingzhi/gomath/stat/correlation"
	"github.com/mingzhi/gomath/stat/desc"
	"math"
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

type MeanVar struct {
	Mean           *desc.Mean
	Var            *desc.Variance
	BiasCorrection bool
}

func NewMeanVar(biasCorrection bool) *MeanVar {
	mv := MeanVar{}
	mv.Mean = desc.NewMean()
	mv.Var = desc.NewVariance()
	mv.BiasCorrection = biasCorrection
	return &mv
}

func (m *MeanVar) Increment(v float64) {
	m.Mean.Increment(v)
	m.Var.Increment(v)
}

func (m *MeanVar) Append(mv2 *MeanVar) {
	m.Mean.Append(mv2.Mean)
	m.Var.Append(mv2.Var)
}

type KsCalculator struct {
	*MeanVar
}

func NewKsCalculator() *KsCalculator {
	kc := KsCalculator{}
	kc.MeanVar = NewMeanVar(true)
	return &kc
}

func (k *KsCalculator) Append(k2 *KsCalculator) {
	k.MeanVar.Append(k2.MeanVar)
}

type MeanCovCalculator struct {
	MeanVars []*MeanVar
}

func NewMeanCovCalculator(maxl int, biasCorrection bool) *MeanCovCalculator {
	s := MeanCovCalculator{}
	s.MeanVars = make([]*MeanVar, maxl)
	for i := 0; i < maxl; i++ {
		s.MeanVars[i] = NewMeanVar(biasCorrection)
	}

	return &s
}

func (s *MeanCovCalculator) Increment(xs, ys []float64, i int) {
	cov := correlation.NewBivariateCovariance(false)
	for i := 0; i < len(xs); i++ {
		x, y := xs[i], ys[i]
		cov.Increment(x, y)
	}
	v := cov.GetResult()
	if !math.IsNaN(v) {
		s.MeanVars[i].Increment(v)
	}
}

func (s *MeanCovCalculator) Append(s2 *MeanCovCalculator) {
	for i := 0; i < len(s2.MeanVars); i++ {
		if len(s.MeanVars) > i {
			s.MeanVars[i].Append(s2.MeanVars[i])
		} else {
			s.MeanVars = append(s.MeanVars, s2.MeanVars[i])
		}
	}
}

type Calculators struct {
	Ks               *KsCalculator
	TCov             *CovCalculator
	SCov, MCov, RCov *MeanCovCalculator
}

func (c *Calculators) Append(c2 *Calculators) {
	c.Ks.Append(c2.Ks)
	c.TCov.Append(c2.TCov)
	c.SCov.Append(c2.SCov)
	c.RCov.Append(c2.RCov)
	c.MCov.Append(c2.MCov)
}

func NewCalculators(maxl int, biasCorrection bool) *Calculators {
	c := Calculators{}
	c.Ks = NewKsCalculator()
	c.TCov = NewCovCalculator(maxl, biasCorrection)
	c.SCov = NewMeanCovCalculator(maxl, biasCorrection)
	c.MCov = NewMeanCovCalculator(maxl, biasCorrection)
	c.RCov = NewMeanCovCalculator(maxl, biasCorrection)
	return &c
}
