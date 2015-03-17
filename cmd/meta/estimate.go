package main

import (
	"fmt"
	"github.com/gonum/plot"
	"github.com/gonum/plot/plotter"
	"github.com/gonum/plot/vg"
	"github.com/mingzhi/gomath/stat/regression"
	"github.com/mingzhi/meta"
	"path/filepath"
	"sort"
	"strings"
)

type Estimation struct {
	FuncName string
	Pos      int
	Strain   meta.Strain
	Ks       float64
	RM       float64
	RSquare  float64
}

// Command to estimate population parameters
// from sequence correlations.
type cmdEstimate struct {
	cmdConfig // embed config parser.
}

// Run command.
func (cmd *cmdEstimate) Run(args []string) {
	cmd.ParseConfig()
	cmd.LoadSpeciesMap()
	MakeDir(filepath.Join(*cmd.workspace, cmd.plotBase))

	INFO.Printf("Fitting Range: %d - %d\n", cmd.fitStart, cmd.fitEnd)
	INFO.Printf("Fitting RSquare Cutoff: %.3f\n", cmd.fitRSquare)

	for prefix, strains := range cmd.speciesMap {
		estimations := []Estimation{}
		for _, funcName := range cmd.covReadsFunctions {
			for _, pos := range cmd.positions {
				if pos != 4 {
					continue
				}
				for _, s := range strains {
					covFileName := fmt.Sprintf("%s_%s_pos%d.json", s.Path, funcName, pos)
					covFilePath := filepath.Join(*cmd.workspace, cmd.covOutBase, covFileName)
					cr := readCovResult(covFilePath)
					simple := regression.NewSimple()

					for i := 0; i < len(cr.CtIndices); i++ {
						index := cr.CtIndices[i]
						if index <= cmd.fitEnd && index >= cmd.fitStart {
							x := float64(index)
							y := 1 / cr.Ct[i]
							simple.Add(x, y)
						}
					}
					b1 := simple.Slope()
					b0 := simple.Intercept()
					ks := cr.Ks
					ratio := b1 * (1.0 + ks*4.0/3.0) / (2.0 * ks * b0)
					est := Estimation{}
					est.Strain = s
					est.FuncName = funcName
					est.Pos = pos
					est.Ks = ks
					est.RM = ratio
					est.RSquare = simple.RSquare()
					estimations = append(estimations, est)
				}
			}
		}
		cmd.Plot(prefix, estimations)
	}
}

type keyType struct {
	FuncName string
	Pos      int
}

type keyTypes []keyType

func (k keyTypes) Len() int {
	return len(k)
}

func (k keyTypes) Swap(i, j int) {
	k[i], k[j] = k[j], k[i]
}

type byFuncName struct{ keyTypes }

func (b byFuncName) Less(i, j int) bool {
	return b.keyTypes[i].FuncName < b.keyTypes[j].FuncName
}

func (cmd *cmdEstimate) Plot(title string, estimations []Estimation) {

	m := make(map[keyType][]float64)
	for _, est := range estimations {
		tp := keyType{est.FuncName, est.Pos}
		if est.RSquare >= cmd.fitRSquare {
			m[tp] = append(m[tp], est.RM)
		}
	}

	// Sort keys by function name.
	kts := keyTypes{}
	for tp, _ := range m {
		kts = append(kts, tp)
	}
	sort.Sort(byFuncName{kts})

	p, err := plot.New()
	if err != nil {
		ERROR.Fatalf("Cannot create plot: %v\n", err)
	}

	p.Title.Text = strings.Replace(title, "_", " ", -1)
	p.Y.Label.Text = "r/m"

	w := vg.Points(20)
	index := 0.0
	names := []string{}
	for _, tp := range kts {
		values := m[tp]
		funcName := strings.Replace(tp.FuncName, "Cov_", "", -1)
		funcName = strings.Replace(funcName, "_", " ", -1)
		name := fmt.Sprintf("%s\nn = %d", funcName, len(values))
		b, err := plotter.NewBoxPlot(w, index, plotter.Values(values))
		if err != nil {
			ERROR.Fatalf("Cannot create boxplot: %v\n", err)
		}
		index++
		names = append(names, name)
		p.Add(b)
	}

	p.NominalX(names...)

	fileName := title + ".pdf"
	filePath := filepath.Join(*cmd.workspace, cmd.plotBase, fileName)
	if err := p.Save(6, 4, filePath); err != nil {
		panic(err)
	}
}
