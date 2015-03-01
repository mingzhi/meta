package meta

import (
	"math"
	"strconv"
	"strings"
)

// Container for blast 6 out result.
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

// Parse and return a Hit struct.
func parseHit(fields []string) Hit {
	h := Hit{}
	h.QSeqid = strings.Split(fields[0], "|")[1]
	h.SSeqid = strings.Split(fields[1], "|")[1]
	h.PIdent = atof(fields[2])
	h.Length = atoi(fields[3])
	h.Mismatch = atoi(fields[4])
	h.GapOpen = atoi(fields[5])
	h.QStart = atoi(fields[6])
	h.QEnd = atoi(fields[7])
	h.SStart = atoi(fields[8])
	h.SEnd = atoi(fields[9])
	h.EValue = atof(fields[10])
	h.BitScore = atof(fields[11])

	return h
}

// String to float64 helper.
func atof(s string) float64 {
	if s == "*" {
		return math.NaN()
	}

	f, err := strconv.ParseFloat(strings.TrimSpace(s), 64)
	if err != nil {
		panic(err)
	}
	return f
}

// String to int helper.
func atoi(s string) int {
	if s == "*" {
		return 0
	}

	i, err := strconv.Atoi(strings.TrimSpace(s))
	if err != nil {
		panic(err)
	}
	return i
}
