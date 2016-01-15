package main

import (
	"github.com/mingzhi/gomath/stat/desc/meanvar"
)

// Covariance contains cov structure
type Covariance struct {
	XY, X, Y float64
	N        int
}

// NewCovariance create a Covariance
func NewCovariance() *Covariance {
	return &Covariance{}
}

// Increment add data to the calculator
func (c *Covariance) Increment(x, y float64) {
	c.XY += x * y
	c.X += x
	c.Y += y
	c.N++
}

// Append merges another covariance.
func (c *Covariance) Append(c1 *Covariance) {
	c.XY += c1.XY
	c.X += c1.X
	c.Y += c1.Y
	c.N += c1.N
}

// GetResult returns the result.
func (c *Covariance) GetResult() float64 {
	var v float64
	v = c.XY/float64(c.N) - (c.X/float64(c.N))*(c.Y/float64(c.N))
	return v
}

// GetN returns N.
func (c *Covariance) GetN() int {
	return c.N
}

// GetMeanX return mean of X.
func (c *Covariance) GetMeanX() float64 {
	return c.X / float64(c.N)
}

// GetMeanY returns mean of Y.
func (c *Covariance) GetMeanY() float64 {
	return c.Y / float64(c.N)
}

// Covariances contains an array of Coveriance.
type Covariances struct {
	corrs []*Covariance
}

// NewCovariances create a new Covariances
func NewCovariances(maxl int) *Covariances {
	cc := Covariances{}
	for i := 0; i < maxl; i++ {
		bc := NewCovariance()
		cc.corrs = append(cc.corrs, bc)
	}
	return &cc
}

// Increment add data (x, y) to the l Covariance.
func (c *Covariances) Increment(l int, x, y float64) {
	c.corrs[l].Increment(x, y)
}

// Append append a Covariance to the l Covariance.
func (c *Covariances) Append(l int, c1 *Covariance) {
	c.corrs[l].Append(c1)
}

// GetResult returns the result.
func (c *Covariances) GetResult(l int) float64 {
	return c.corrs[l].GetResult()
}

// GetN returns the number of data points.
func (c *Covariances) GetN(l int) int {
	return c.corrs[l].GetN()
}

// GetMeanX returns x_bar.
func (c *Covariances) GetMeanX(l int) float64 {
	return c.corrs[l].GetMeanX()
}

// GetMeanY returns y_bar.
func (c *Covariances) GetMeanY(l int) float64 {
	return c.corrs[l].GetMeanY()
}

// MeanVariances is an array of MeanVar.
type MeanVariances struct {
	meanvars []*meanvar.MeanVar
}

// NewMeanVariances return a new MeanVariances of length maxl.
func NewMeanVariances(size int) *MeanVariances {
	mcc := MeanVariances{}
	for i := 0; i < size; i++ {
		mv := meanvar.New()
		mcc.meanvars = append(mcc.meanvars, mv)
	}
	return &mcc
}

// Increment add a data point to the lst MeanVar.
func (m *MeanVariances) Increment(l int, v float64) {
	m.meanvars[l].Increment(v)
}

// GetMean returns the mean of the lst MeanVar.
func (m *MeanVariances) GetMean(l int) float64 {
	return m.meanvars[l].Mean.GetResult()
}

// GetVar returns the variance of the lst MeanVar.
func (m *MeanVariances) GetVar(l int) float64 {
	return m.meanvars[l].Var.GetResult()
}

// GetN returns the size of the data points of the lst MeanVar.
func (m *MeanVariances) GetN(l int) int {
	return m.meanvars[l].Mean.GetN()
}

// Calculator contains individual calculators.
type Calculator struct {
	MaxL int
	Cs   *MeanVariances
	Cr   *Covariances
	Ct   *Covariances
}

// NewCalculator returns a new Calculator
func NewCalculator(maxl int) *Calculator {
	c := Calculator{}
	c.MaxL = maxl
	c.Cs = NewMeanVariances(maxl)
	c.Cr = NewCovariances(maxl)
	c.Ct = NewCovariances(maxl)
	return &c
}

// Calc add two SNP.
func (c *Calculator) Calc(s1, s2 *SNP) {
	if s1.Pos > s2.Pos {
		s1, s2 = s2, s1
	}

	l := s2.Pos - s1.Pos
	if l < c.MaxL {
		basePairs := pairBases(s1, s2)
		if len(basePairs) >= minDepth {
			cov := covPairedBases(basePairs)
			c.Cs.Increment(l, cov.GetResult())
			c.Ct.Append(l, cov)
			c.Cr.Increment(l, cov.GetMeanX(), cov.GetMeanY())
		}
	}
}

func covPairedBases(pairs []basePair) (c *Covariance) {
	c = NewCovariance()
	for i := 0; i < len(pairs); i++ {
		p1 := pairs[i]
		for j := i + 1; j < len(pairs); j++ {
			p2 := pairs[j]
			var x, y float64
			if p1.a.Base != p2.a.Base {
				x = 1.0
			}

			if p1.b.Base != p2.b.Base {
				y = 1.0
			}

			c.Increment(x, y)
		}
	}
	return c
}

// basePair is a pair of bases, which come from the same read (or paired-end).
type basePair struct {
	a, b *Base
}

func pairBases(s1, s2 *SNP) (pairs []basePair) {
	m := make(map[string]*Base)
	for i := 0; i < len(s1.Bases); i++ {
		b := s1.Bases[i]
		m[b.ReadId] = b
	}

	for i := 0; i < len(s2.Bases); i++ {
		b2 := s2.Bases[i]
		b1, found := m[b2.ReadId]
		if found {
			if b1.Pos > b2.Pos {
				b1, b2 = b2, b1
			}
			pairs = append(pairs, basePair{b1, b2})
		}
	}

	return
}
