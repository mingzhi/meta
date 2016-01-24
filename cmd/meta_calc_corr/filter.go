package main

import (
	"github.com/biogo/hts/sam"
	"log"
	"strings"
)

func FilterReads(input chan *sam.Record) (output chan *sam.Record) {
	output = make(chan *sam.Record)
	go func() {
		defer close(output)
		m := make(map[string]*sam.Record)
		deletes := make(map[string]bool)
		for r := range input {
			if isProperPair(r) {
				if isOnlyMatched(r) && int(r.MapQ) >= minMQ {
					mate, found := m[r.Name]
					if found {
						if isOverlapMatched(r, mate) {
							output <- r
							output <- mate
						}
						delete(m, r.Name)
					} else {
						if deletes[r.Name] {
							delete(deletes, r.Name)
						} else {
							m[r.Name] = r
						}
					}
				} else {
					if deletes[r.Name] {
						delete(deletes, r.Name)
					} else {
						if _, found := m[r.Name]; found {
							delete(m, r.Name)
						} else {
							deletes[r.Name] = true
						}
					}
				}
			}
		}
		log.Println("Finished filtering reads!")
	}()

	return
}

func isProperPair(r *sam.Record) bool {
	return strings.Contains(r.Flags.String(), "P")
}

func isOnlyMatched(r *sam.Record) bool {
	return len(r.Cigar) == 1 && r.Cigar[0].Type() == sam.CigarMatch
}

func isOverlapMatched(r, mate *sam.Record) bool {
	if r.Pos > mate.Pos {
		r, mate = mate, r
	}
	s1, _ := Map2Ref(r)
	s2, _ := Map2Ref(mate)
	matched := true
	for i := 0; i < len(s2); i++ {
		j := mate.Pos - r.Pos + i
		if j >= len(s1) {
			break
		}

		if s1[j] != s2[i] {
			matched = false
			break
		}
	}

	return matched
}
