package main

import (
	"log"
	"os"
)

func MakeDir(d string) {
	err := os.MkdirAll(d, 0666)
	if err != nil {
		log.Panic(err)
	}
}
