package meta

import (
	"strconv"
	"strings"
)

type Hit struct {
	QSeqid   string
	SSeqid   string
	PIdent   float64
	Length   int
	Mismatch int
	GapOpen  int
	QLen     int
	QStart   int
	QEnd     int
	SLen     int
	SStart   int
	SEnd     int
	EValue   float64
	BitScore float64
}

func parseHit(fields []string) Hit {
	h := Hit{}
	h.QSeqid = strings.Split(fields[0], "|")[1]
	h.SSeqid = strings.Split(fields[1], "|")[1]
	h.PIdent = atof(fields[2])
	h.Length = atoi(fields[3])
	h.Mismatch = atoi(fields[4])
	h.GapOpen = atoi(fields[5])
	h.QLen = atoi(fields[6])
	h.QStart = atoi(fields[7])
	h.QEnd = atoi(fields[8])
	h.SLen = atoi(fields[9])
	h.SStart = atoi(fields[10])
	h.SEnd = atoi(fields[11])
	h.EValue = atof(fields[12])
	h.BitScore = atof(fields[13])
	return h
}

func atof(s string) float64 {
	f, err := strconv.ParseFloat(strings.TrimSpace(s), 64)
	if err != nil {
		panic(err)
	}
	return f
}

func atoi(s string) int {
	i, err := strconv.Atoi(strings.TrimSpace(s))
	if err != nil {
		panic(err)
	}
	return i
}
