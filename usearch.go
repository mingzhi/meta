package meta

import (
	"bytes"
	"encoding/csv"
	"io/ioutil"
	"log"
	"os"
	"os/exec"
)

func UsearchMakeUDB(f string) {
	// Info.Println(f)
	cmd := exec.Command("usearch", "-makeudb_usearch", f, "-output", f+".udb")
	run(cmd)

	return
}

func UsearchGlobal(q, db string) []Hit {
	// Prepare temp file for output.
	temp, err := ioutil.TempFile("", "usearchglobal")
	if err != nil {
		log.Panic(err)
	}
	temp.Close()

	// Prepare command and run it.
	cmd := exec.Command("usearch", "-usearch_global", q, "-db", db, "-id", "0.8", "-blast6out", temp.Name(), "-top_hit_only")
	run(cmd)

	// Open and read results.
	f, err := os.Open(temp.Name())
	if err != nil {
		log.Panic(err)
	}
	defer f.Close()

	csvRd := csv.NewReader(f)
	csvRd.Comma = '\t'
	rows, err := csvRd.ReadAll()
	if err != nil {
		log.Panic(err)
	}

	// Parse hits.
	hits := []Hit{}
	for _, fields := range rows {
		h := parseHit(fields)
		hits = append(hits, h)
	}

	return hits
}

func run(cmd *exec.Cmd) {
	stderr := new(bytes.Buffer)
	cmd.Stderr = stderr
	if err := cmd.Run(); err != nil {
		log.Fatalf("Error when %s: %s\n", cmd.Args[1], stderr.String())
	}

	return
}

func IsUsearchDBExist(f string) bool {
	fileName := f + ".udb"
	if _, err := os.Stat(fileName); err == nil {
		return true
	} else {
		return false
	}
}
