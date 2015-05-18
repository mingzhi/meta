package main

import (
	"encoding/json"
	"fmt"
	"github.com/mingzhi/meta/fit"
	"github.com/mingzhi/meta/strain"
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
					results := fromJson(filePath)
					fitResults := fitExp(results, cmd.fitStart, cmd.fitEnd)
					fitFileOutPath := filepath.Join(*cmd.workspace, cmd.fitOutBase, s.Path, filePrefix+"_boot.json")
					toJson(fitFileOutPath, fitResults)
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

func fitExp(results []CovResult, fitStart, fitEnd int) (fitResults []FitResult) {
	jobs := make(chan CovResult)
	go func() {
		defer close(jobs)
		for _, r := range results {
			jobs <- r
		}
	}()

	ncpu := runtime.GOMAXPROCS(0)
	done := make(chan bool)
	fitResChan := make(chan FitResult)
	for i := 0; i < ncpu; i++ {
		go func() {
			for r := range jobs {
				fr := FitResult{}
				fr.Ks = r.Ks
				xdata := []float64{}
				for i := fitStart; i < fitEnd; i++ {
					xdata = append(xdata, float64(r.CtIndices[i]))
				}
				ydata := r.Ct[fitStart:fitEnd]
				par := fit.FitExp(xdata, ydata)
				fr.B0 = par[0]
				fr.B1 = par[1]
				fr.B2 = par[2]
				fitResChan <- fr
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

	for fr := range fitResChan {
		if !isNaN(fr) {
			fitResults = append(fitResults, fr)
		}
	}

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

func fromJson(filePath string) (results []CovResult) {
	f, err := os.Open(filePath)
	if err != nil {
		panic(err)
	}
	defer f.Close()

	d := json.NewDecoder(f)
	if err := d.Decode(&results); err != nil {
		panic(err)
	}

	return
}

func toJson(filePath string, fitResults []FitResult) {
	f, err := os.Create(filePath)
	if err != nil {
		panic(err)
	}
	defer f.Close()

	e := json.NewEncoder(f)
	if err := e.Encode(fitResults); err != nil {
		panic(err)
	}
}
