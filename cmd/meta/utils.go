package main

import (
	"encoding/json"
	"github.com/mingzhi/meta"
	"log"
	"os"
)

type CovResult struct {
	Ks, VarKs float64
	Ct        []float64
	CtIndices []int
	CtN       []int
	N         int
}

func MakeDir(d string) {
	err := os.MkdirAll(d, 0666)
	if err != nil {
		ERROR.Fatalln(err)
	}
}

func registerLogger() {
	meta.Info = log.New(os.Stdout, "INFO: ", log.Ldate|log.Ltime|log.Lshortfile)
	meta.Warn = log.New(os.Stdout, "WARN: ", log.Ldate|log.Ltime|log.Lshortfile)
}

// Save cov result to a json file.
func save2Json(cr CovResult, fileName string) {
	f, err := os.Create(fileName)
	if err != nil {
		ERROR.Fatalln(err)
	}
	defer f.Close()

	ec := json.NewEncoder(f)
	if err := ec.Encode(cr); err != nil {
		ERROR.Fatalln(err)
	}
}
