package main

import (
	"bytes"
	"github.com/mingzhi/gomath/stat/correlation"
	"github.com/mingzhi/gomath/stat/desc"
	"github.com/mingzhi/gomath/stat/desc/meanvar"
	"math"
)

type CovCalculator struct {
	corrs []*correlation.BivariateCovariance
}

func NewCovCalculator(maxl int, bias bool) *CovCalculator {
	cc := CovCalculator{}
	for i := 0; i < maxl; i++ {
		bc := correlation.NewBivariateCovariance(bias)
		cc.corrs = append(cc.corrs, bc)
	}
	return &cc
}

func (c *CovCalculator) Increment(l int, x, y float64) {
	c.corrs[l].Increment(x, y)
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
}

func NewCalculator(maxl int) *Calculator {
	c := Calculator{}
	c.MaxL = maxl
	bias := false
	c.Cs = NewMeanCovCalculator(maxl)
	c.Cr = NewCovCalculator(maxl, bias)
	c.Ks = NewKsCalculator()
	return &c
}

func (c *Calculator) Calc(s1, s2 *SNP) {
	if s1.Pos > s2.Pos {
		s1, s2 = s2, s1
	}

	l := s2.Pos - s1.Pos
	if l < c.MaxL {
		x, nx := calcPi(s1.Bases)
		y, ny := calcPi(s2.Bases)
		if nx > minDepth && ny > minDepth {
			c.Cr.Increment(l, x, y)
			cov := covSNPs(s1, s2)
			if cov.GetN() > minDepth {
				c.Cs.Increment(l, cov.GetResult())
			}
			c.Ks.Increment(x)
			c.Ks.Increment(y)
		}
	}
}

func calcPi(bases []*Base) (pi float64, n int) {
	bs := []byte{}
	for i := 0; i < len(bases); i++ {
		bs = append(bs, bases[i].Base)
	}
	// convert bases to upper case.
	upperBases := bytes.ToUpper(bs)

	mean := desc.NewMean()
	for i := 0; i < len(upperBases); i++ {
		if upperBases[i] == '*' {
			continue
		}

		for j := i + 1; j < len(upperBases); j++ {
			if upperBases[j] == '*' {
				continue
			}

			if upperBases[i] != upperBases[j] {
				mean.Increment(1.0)
			} else {
				mean.Increment(0.0)
			}
		}
	}

	pi = mean.GetResult()
	n = mean.GetN()
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

func covSNPs(s1, s2 *SNP) (c *correlation.BivariateCovariance) {
	c = correlation.NewBivariateCovariance(false)
	pairs := pairBases(s1, s2)
	for i := 0; i < len(pairs); i++ {
		p1 := pairs[i]
		for j := i + 1; j < len(pairs); j++ {
			p2 := pairs[j]
			var x, y float64
			if p1.A.Base != p2.A.Base {
				x = 1.0
			} else {
				x = 0.0
			}

			if p1.B.Base != p2.B.Base {
				y = 1.0
			} else {
				y = 0.0
			}
			c.Increment(x, y)
		}
	}

	return c
}

func Calc(snpChan chan *SNP, profile []byte, t byte, maxl, geneLength int) (cChan chan *Calculator) {

	storage := []*SNP{}

	cChan = make(chan *Calculator)

	go func() {
		defer close(cChan)

		totalLength := geneLength
		var c *Calculator
		c = NewCalculator(maxl)
		for snp := range snpChan {
			if snp.Pos >= totalLength {
				cChan <- c
				c = NewCalculator(maxl)
				totalLength += geneLength
			}

			storage = append(storage, snp)

			if len(storage) > maxl {
				s1 := storage[0]
				if checkPosType(t, profile[s1.Pos]) {
					for i := 1; i < len(storage); i++ {
						s2 := storage[i]
						if checkPosType(t, profile[s2.Pos]) {
							c.Calc(s1, s2)
						}
					}
				}
				storage = storage[1:]
			}
		}

		for i := 0; i < len(storage); i++ {
			s1 := storage[i]
			if profile[s1.Pos] == t {
				for j := i + 1; j < len(storage); j++ {
					s2 := storage[j]
					if profile[s2.Pos] == t {
						c.Calc(s1, s2)
					}
				}
			}

		}

		cChan <- c
	}()

	return
}

func Collect(maxl int, cChan chan *Calculator) (means, covs, ks []*meanvar.MeanVar) {
	for i := 0; i < maxl; i++ {
		means = append(means, meanvar.New())
		covs = append(covs, meanvar.New())
	}

	ks = make([]*meanvar.MeanVar, 2)
	ks[0] = meanvar.New()
	ks[1] = meanvar.New()

	for c := range cChan {
		cr := c.Cr
		for i := 0; i < c.MaxL; i++ {
			v := cr.GetResult(i)
			if !math.IsNaN(v) {
				covs[i].Increment(v)
			}
		}

		cs := c.Cs
		for i := 0; i < c.MaxL; i++ {
			v := cs.GetMean(i)
			if !math.IsNaN(v) {
				means[i].Increment(v)
			}
		}

		ksm := c.Ks.GetMean()
		ksv := c.Ks.GetVariance()
		if !math.IsNaN(ksv) {
			ks[0].Increment(ksm)
			ks[1].Increment(ksv)
		}
	}

	return
}

func checkPosType(t, t1 byte) bool {
	if t == ThirdPos {
		if t1 == ThirdPos || t1 == FourFold {
			return true
		} else {
			return false
		}
	} else {
		return t == t1
	}
}
