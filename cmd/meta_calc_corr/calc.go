package main

import (
	"fmt"
	"github.com/mingzhi/gomath/stat/desc/meanvar"
	"github.com/mingzhi/ncbiftp/genomes/profiling"
	"math"
	"runtime"
)

// Calc perform calculations.
func Calc(snpChan chan *SNP, profile []profiling.Pos, posType byte, maxl int) (cChan chan *Calculator) {
	ncpu := runtime.GOMAXPROCS(0)
	// create job chan
	// each job is a list of SNP in a gene.
	geneSNPChan := make(chan []*SNP, ncpu)
	go func() {
		defer close(geneSNPChan)
		var geneName string
		var storage []*SNP
		for snp := range snpChan {
			if checkPosType(posType, profile[snp.Pos-1].Type) {
				geneName1 := profile[snp.Pos-1].Gene
				if len(storage) != 0 && geneName != geneName1 {
					geneSNPChan <- storage
					storage = []*SNP{}
					geneName = geneName1
					fmt.Println(geneName)
				}
				storage = append(storage, snp)
			}
		}

		if len(storage) != 0 {
			geneSNPChan <- storage
		}
	}()

	cChan = make(chan *Calculator)
	done := make(chan bool)
	for i := 0; i < ncpu; i++ {
		go func() {
			for storage := range geneSNPChan {
				c := NewCalculator(maxl)
				calc(c, storage, posType, profile, maxl)
				cChan <- c
			}
			done <- true
		}()
	}

	go func() {
		defer close(cChan)
		for i := 0; i < ncpu; i++ {
			<-done
		}
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

func Collect(maxl int, cChan chan *Calculator) (means, covs, xbars, ybars, totals []*meanvar.MeanVar) {
	for i := 0; i < maxl; i++ {
		means = append(means, meanvar.New())
		covs = append(covs, meanvar.New())
		totals = append(totals, meanvar.New())
		xbars = append(xbars, meanvar.New())
		ybars = append(ybars, meanvar.New())
	}

	for c := range cChan {

		for i := 0; i < c.MaxL; i++ {
			cr := c.Cr.GetResult(i)
			cs := c.Cs.GetMean(i)
			ct := c.Ct.GetResult(i)
			xbar := c.Cr.GetMeanX(i)
			ybar := c.Cr.GetMeanY(i)
			n := c.Cr.GetN(i)

			if !math.IsNaN(cr) && !math.IsNaN(cs) && !math.IsNaN(ct) && n >= 50 {
				covs[i].Increment(cr)
				means[i].Increment(cs)
				totals[i].Increment(ct)
				xbars[i].Increment(xbar)
				ybars[i].Increment(ybar)
			}
		}
	}

	return
}

func checkPosType(t, t1 byte) bool {
	isFirstPos := t1 == profiling.FirstPos
	isSecondPos := t1 == profiling.SecondPos
	isThirdPos := t1 == profiling.ThirdPos
	isFourFold := t1 == profiling.FourFold

	if t == profiling.Coding {
		if isFirstPos || isSecondPos || isThirdPos || isFourFold {
			return true
		}
		return false
	}

	if t == profiling.ThirdPos {
		if isThirdPos || isFourFold {
			return true
		}
		return false
	}

	return t == t1
}
