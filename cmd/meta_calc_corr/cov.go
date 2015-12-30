package main

import (
	"bytes"
	"github.com/mingzhi/gomath/stat/desc"
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

type CovCalculator struct {
	corrs []*Covariance
}

func NewCovCalculator(maxl int) *CovCalculator {
	cc := CovCalculator{}
	for i := 0; i < maxl; i++ {
		bc := NewCovariance()
		cc.corrs = append(cc.corrs, bc)
	}
	return &cc
}

func (c *CovCalculator) Increment(l int, x, y float64) {
	c.corrs[l].Increment(x, y)
}

func (c *CovCalculator) Append(l int, c1 *Covariance) {
	c.corrs[l].Append(c1)
}

func (c *CovCalculator) GetResult(l int) float64 {
	return c.corrs[l].GetResult()
}

func (c *CovCalculator) GetN(l int) int {
	return c.corrs[l].GetN()
}

type MeanCovCalculator struct {
	meanvars []*meanvar.MeanVar
}

func NewMeanCovCalculator(maxl int) *MeanCovCalculator {
	mcc := MeanCovCalculator{}
	for i := 0; i < maxl; i++ {
		mv := meanvar.New()
		mcc.meanvars = append(mcc.meanvars, mv)
	}
	return &mcc
}

func (m *MeanCovCalculator) Increment(l int, v float64) {
	m.meanvars[l].Increment(v)
}

func (m *MeanCovCalculator) GetMean(l int) float64 {
	return m.meanvars[l].Mean.GetResult()
}

func (m *MeanCovCalculator) GetVar(l int) float64 {
	return m.meanvars[l].Var.GetResult()
}

func (m *MeanCovCalculator) GetN(l int) int {
	return m.meanvars[l].Mean.GetN()
}

type KsCalculator struct {
	mean     *desc.Mean
	variance *desc.Variance
}

func NewKsCalculator() *KsCalculator {
	ks := KsCalculator{}
	ks.mean = desc.NewMean()
	ks.variance = desc.NewVariance()
	return &ks
}

func (k *KsCalculator) Increment(f float64) {
	k.mean.Increment(f)
	k.variance.Increment(f)
}

func (k *KsCalculator) GetMean() float64 {
	return k.mean.GetResult()
}

func (k *KsCalculator) GetVariance() float64 {
	return k.variance.GetResult()
}

func (k *KsCalculator) GetN() int {
	return k.mean.GetN()
}

type Calculator struct {
	MaxL int
	Cs   *MeanCovCalculator
	Cr   *CovCalculator
	Ks   *KsCalculator
	Ct   *CovCalculator
}

func NewCalculator(maxl int) *Calculator {
	c := Calculator{}
	c.MaxL = maxl
	c.Cs = NewMeanCovCalculator(maxl)
	c.Cr = NewCovCalculator(maxl)
	c.Ks = NewKsCalculator()
	c.Ct = NewCovCalculator(maxl)
	return &c
}

func (c *Calculator) Calc(s1, s2 *SNP) {
	if s1.Pos > s2.Pos {
		s1, s2 = s2, s1
	}

	l := s2.Pos - s1.Pos
	if l < c.MaxL {
		cov := covSNPs(s1, s2)
		cutoff := minDepth * (minDepth - 1) / 2
		if cov.GetN() >= cutoff {
			c.Cs.Increment(l, cov.GetResult())
			c.Ct.Append(l, cov)
			c.Cr.Increment(l, cov.GetMeanX(), cov.GetMeanY())
			c.Ks.Increment(cov.GetMeanX())
			c.Ks.Increment(cov.GetMeanY())
		}
	}
}

func calcPi(bases []*Base) (pi float64, depth int) {
	bs := []byte{}
	for i := 0; i < len(bases); i++ {
		bs = append(bs, bases[i].Base)
	}

	counts := make([]int, 4)
	ss := bytes.ToUpper(bs)
	for i := 0; i < len(ss); i++ {
		b := ss[i]
		switch b {
		case 'A':
			counts[0]++
			break
		case 'T':
			counts[1]++
			break
		case 'C':
			counts[2]++
			break
		case 'G':
			counts[3]++
		}
	}

	total := 0
	cross := 0
	for i := 0; i < len(counts); i++ {
		x := counts[i]
		for j := i + 1; j < len(counts); j++ {
			y := counts[j]
			cross += x * y
		}
		total += x

	}

	pi = float64(cross) / float64(total*(total-1)/2)
	depth = total

	return
}

type BasePair struct {
	A, B *Base
}

func pairBases(s1, s2 *SNP) (pairs []BasePair) {
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
			pairs = append(pairs, BasePair{b1, b2})
		}
	}

	return
}

func covSNPs(s1, s2 *SNP) (c *Covariance) {
	c = NewCovariance()
	pairs := pairBases(s1, s2)
	for i := 0; i < len(pairs); i++ {
		p1 := pairs[i]
		for j := i + 1; j < len(pairs); j++ {
			p2 := pairs[j]
			var x, y float64
			if p1.A.Base != p2.A.Base {
				x = 1.0
			}

			if p1.B.Base != p2.B.Base {
				y = 1.0
			}

			c.Increment(x, y)
		}
	}

	return c
}
