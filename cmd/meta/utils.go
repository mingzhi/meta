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
	NReads    int
}

func MakeDir(d string) {
	err := os.MkdirAll(d, 0777)
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
		ERROR.Panicln(err)
	}
	defer f.Close()

	ec := json.NewEncoder(f)
	if err := ec.Encode(cr); err != nil {
		ERROR.Println(cr)
		ERROR.Printf("Ks: %f\n", cr.Ks)
		ERROR.Printf("VarKs: %f\n", cr.VarKs)
		ERROR.Panicln(err)
	}
}

// Read cov result.
func readCovResult(fileName string) (cr CovResult) {
	r, err := os.Open(fileName)
	if err != nil {
		ERROR.Fatalln("Cannot open file: %s, %v", fileName, err)
	}
	defer r.Close()

	decoder := json.NewDecoder(r)
	err = decoder.Decode(&cr)
	if err != nil {
		ERROR.Fatalln("Cannot decode file: %s, with CovResult: %v", fileName, err)
	}
	return
}
