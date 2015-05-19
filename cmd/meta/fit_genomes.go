package main

import (
	"encoding/json"
	"fmt"
	"github.com/mingzhi/meta/fit"
	"github.com/mingzhi/meta/strain"
	"io"
	"log"
	"math"
	"os"
	"path/filepath"
	"runtime"
)

type cmdFitGenomes struct {
	cmdConfig
}

func (cmd *cmdFitGenomes) Init() {
	// Parse config and settings.
	cmd.ParseConfig()
	// Load species map.
	cmd.LoadSpeciesMap()
	// Make output directory.
	MakeDir(filepath.Join(*cmd.workspace, cmd.fitOutBase))
	// Check profile positions.
	if len(cmd.positions) == 0 {
		WARN.Println("Use default position: 4!")
		cmd.positions = append(cmd.positions, 4)
	}
	runtime.GOMAXPROCS(runtime.NumCPU())
}

func (cmd *cmdFitGenomes) Run(args []string) {
	cmd.Init()
	type job struct {
		strains []strain.Strain
		pos     int
		typ     string
		funcT   string
	}
	jobs := make(chan job)
	go func() {
		defer close(jobs)
		for _, strains := range cmd.speciesMap {
			for _, pos := range cmd.positions {
				for _, name := range []string{"core", "disp", "pan"} {
					for _, funcType := range []string{"Cov_Genomes_vs_Genome", "Cov_Genomes_vs_Genomes"} {
						j := job{}
						j.strains = strains
						j.pos = pos
						j.typ = name
						j.funcT = funcType
						jobs <- j
					}
				}
			}
		}
	}()

	ncpu := runtime.GOMAXPROCS(0)
	done := make(chan bool)
	for i := 0; i < ncpu; i++ {
		go func() {
			for j := range jobs {
				pos := j.pos
				name := j.typ
				funcType := j.funcT
				strains := j.strains
				cmd.RunOne(strains, pos, name, funcType)
			}
			done <- true
		}()
	}

	for i := 0; i < ncpu; i++ {
		<-done
	}
}

func (cmd *cmdFitGenomes) RunOne(strains []strain.Strain, pos int, name string, funcType string) {
	jobs := make(chan strain.Strain)
	go func() {
		defer close(jobs)
		for _, s := range strains {
			jobs <- s
		}
	}()

	ncpu := runtime.GOMAXPROCS(0)
	done := make(chan bool)
	for i := 0; i < ncpu; i++ {
		go func() {
			for s := range jobs {
				MakeDir(filepath.Join(*cmd.workspace, cmd.fitOutBase, s.Path))
				for _, g := range s.Genomes {
					filePrefix := fmt.Sprintf("%s_%s_%s_pos%d", g.RefAcc(), funcType, name, pos)
					filePath := filepath.Join(*cmd.workspace, cmd.covOutBase, s.Path, filePrefix+"_boot.json")
					var f fitFunc
					for _, fitCon := range cmd.fitControls {
						if fitCon.end-fitCon.start > 0 {
							name := fitCon.name
							switch name {
							case "exp":
								f = fitExp
							case "hyper":
								f = fitHyper
							default:
								f = nil
							}

							if f != nil {
								resChan := fromJson(filePath)
								fitResChan := doFit(f, resChan, fitCon.start, fitCon.end)
								fitFileOutPath := filepath.Join(*cmd.workspace, cmd.fitOutBase, s.Path, filePrefix+"_"+name+"_boot.json")
								toJson(fitFileOutPath, fitResChan)
							}
						}
					}
				}
			}
			done <- true
		}()
	}

	for i := 0; i < ncpu; i++ {
		<-done
	}
}

type FitResult struct {
	Ks         float64
	B0, B1, B2 float64
}

type fitFunc func(xdata, ydata []float64) FitResult

func doFit(f fitFunc, resChan chan CovResult, fitStart, fitEnd int) (fitResChan chan FitResult) {
	ncpu := runtime.GOMAXPROCS(0)
	done := make(chan bool)
	fitResChan = make(chan FitResult)
	for i := 0; i < ncpu; i++ {
		go func() {
			for r := range resChan {
				fr := FitResult{}
				fr.Ks = r.Ks
				xdata := []float64{}
				ydata := []float64{}
				for i := 0; i < len(r.CtIndices) && r.CtIndices[i] < fitEnd; i++ {
					if r.CtIndices[i] >= fitStart {
						xdata = append(xdata, float64(r.CtIndices[i]))
						ydata = append(ydata, r.Ct[i])
					}
				}

				fitResChan <- f(xdata, ydata)
			}
			done <- true
		}()
	}

	go func() {
		defer close(fitResChan)
		for i := 0; i < ncpu; i++ {
			<-done
		}
	}()

	return
}

func fitHyper(xdata, ydata []float64) (res FitResult) {
	par := fit.FitHyper(xdata, ydata)
	res.B0 = par[0]
	res.B1 = par[1]
	return
}

func fitExp(xdata, ydata []float64) (res FitResult) {
	par := fit.FitExp(xdata, ydata)
	res.B0 = par[0]
	res.B1 = par[1]
	res.B2 = par[2]
	return
}

func isNaN(res FitResult) bool {
	floats := []float64{}
	floats = append(floats, res.Ks)
	floats = append(floats, res.B1)
	floats = append(floats, res.B0)
	floats = append(floats, res.B2)

	for _, v := range floats {
		if math.IsNaN(v) {
			return true
		}
	}

	return false
}

func fromJson(filePath string) (resChan chan CovResult) {
	f, err := os.Open(filePath)
	if err != nil {
		panic(err)
	}

	d := json.NewDecoder(f)

	resChan = make(chan CovResult)
	go func() {
		defer close(resChan)
		for {
			res := CovResult{}
			if err := d.Decode(&res); err != nil {
				if err != io.EOF {
					log.Panicln(err)
				}
				break
			} else {
				resChan <- res
			}
		}
		f.Close()
	}()

	return
}

func toJson(filePath string, fitResChan chan FitResult) {
	f, err := os.Create(filePath)
	if err != nil {
		panic(err)
	}
	defer f.Close()

	fitResults := []FitResult{}
	for res := range fitResChan {
		fitResults = append(fitResults, res)
	}

	e := json.NewEncoder(f)
	if err := e.Encode(fitResults); err != nil {
		panic(err)
	}
}
