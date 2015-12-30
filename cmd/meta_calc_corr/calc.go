package main

import (
	"github.com/mingzhi/gomath/stat/desc/meanvar"
	"github.com/mingzhi/ncbiftp/genomes/profiling"
	"math"
)

func Calc(snpChan chan *SNP, profile []profiling.Pos, posType byte, maxl int) (cChan chan *Calculator) {

	cChan = make(chan *Calculator)

	go func() {
		defer close(cChan)

		var storage []*SNP
		var geneName string
		var c *Calculator
		c = NewCalculator(maxl)
		for snp := range snpChan {
			if checkPosType(posType, profile[snp.Pos-1].Type) {
				geneName1 := profile[snp.Pos-1].Gene
				if len(storage) != 0 && geneName != geneName1 {
					calc(c, storage, posType, profile, maxl)
					cChan <- c
					c = NewCalculator(maxl)
					geneName = geneName1
					storage = []*SNP{}

				}
				storage = append(storage, snp)
			}
		}

		calc(c, storage, posType, profile, maxl)

		cChan <- c
	}()

	return
}

func calc(c *Calculator, snpArr []*SNP, t byte, profile []profiling.Pos, maxl int) {
	for i := 0; i < len(snpArr); i++ {
		s1 := snpArr[i]
		if checkPosType(t, profile[s1.Pos-1].Type) {
			for j := i; j < len(snpArr); j++ {
				s2 := snpArr[j]
				if checkPosType(t, profile[s2.Pos-1].Type) {
					c.Calc(s1, s2)
				}
			}
		}
	}
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
			n := cr.GetN(i)
			if !math.IsNaN(v) && n >= 50 {
				covs[i].Increment(v)
			}
		}

		cs := c.Cs
		for i := 0; i < c.MaxL; i++ {
			v := cs.GetMean(i)
			n := cs.GetN(i)
			if !math.IsNaN(v) && n >= 50 {
				means[i].Increment(v)
			}
		}

		ksm := c.Ks.GetMean()
		ksv := c.Ks.GetVariance()
		n := c.Ks.GetN()
		if !math.IsNaN(ksv) && n >= 50 {
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
		}
		return false
	}

	return t == t1
}
