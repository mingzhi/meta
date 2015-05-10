package fit

/*
#include "fit.h"
*/
import "C"
import (
	"fmt"
	"github.com/mingzhi/gomath/stat/regression"

	"unsafe"
)

func FitHyper(t, y []float64) []float64 {
	s := regression.NewSimple()
	for i := 0; i < 10; i++ {
		if i >= len(y) {
			break
		}
		iy := 1.0 / y[i]
		s.Add(t[i], iy)
	}

	n := 2
	m := len(t)
	par := []float64{s.Intercept(), s.Slope()}
	C.fitHyper(C.int(n), (*C.double)(unsafe.Pointer(&par[0])), C.int(m), (*C.double)(unsafe.Pointer(&t[0])), (*C.double)(unsafe.Pointer(&y[0])))

	return par
}

func FitExp(t, y []float64) []float64 {
	n := 3
	m := len(t)
	l := 6
	if len(t) < l {
		l = len(t)
	}
	par := FitHyper(t[:l], y[:l])
	par[1] = par[1] * 100
	par = append(par, 100.0)
	fmt.Println(par)
	C.fitExp(C.int(n), (*C.double)(unsafe.Pointer(&par[0])), C.int(m), (*C.double)(unsafe.Pointer(&t[0])), (*C.double)(unsafe.Pointer(&y[0])))

	return par
}
