package main

import (
	"github.com/mingzhi/meta"
	"log"
	"os"
)

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
