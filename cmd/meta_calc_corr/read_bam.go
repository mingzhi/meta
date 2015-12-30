package main

import (
	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
	"io"
	"os"
)

// ReadBamFile reads bam file, and return the header and a channel of sam records.
func ReadBamFile(fileName string) (h *sam.Header, c chan *sam.Record) {
	// Initialize the channel of sam records.
	c = make(chan *sam.Record)

	// Create a new go routine to read the records.
	go func() {
		// Close the record channel when finished.
		defer close(c)

		// Open file stream, and close it when finished.
		f, err := os.Open(fileName)
		if err != nil {
			panic(err)
		}
		defer f.Close()

		// Create a bam reader, and close it when finished.
		rd := 0 // If rd is zero concurrency is set to GOMAXPROCS.
		bamReader, err := bam.NewReader(f, rd)
		if err != nil {
			panic(err)
		}
		defer bamReader.Close()

		// Read and assign header.
		h = bamReader.Header()

		// Read sam records and send them to the channel,
		// until it hit an error, which raises a panic
		// if it is not a IO EOF.
		for {
			rec, err := bamReader.Read()
			if err != nil {
				if err != io.EOF {
					panic(err)
				}
				break
			}
			c <- rec
		}
	}()

	return
}

// Map2Ref Obtains a read mapping to the reference genome.
func Map2Ref(r *sam.Record) (s []byte, q []byte) {
	p := 0                 // position in the read sequence.
	read := r.Seq.Expand() // read sequence.
	qual := r.Qual
	for _, c := range r.Cigar {
		switch c.Type() {
		case sam.CigarMatch, sam.CigarMismatch, sam.CigarEqual:
			s = append(s, read[p:p+c.Len()]...)
			q = append(q, qual[p:p+c.Len()]...)
			p += c.Len()
		case sam.CigarInsertion, sam.CigarSoftClipped, sam.CigarHardClipped:
			p += c.Len()

		case sam.CigarDeletion, sam.CigarSkipped:
			for i := 0; i < c.Len(); i++ {
				s = append(s, '*')
				q = append(q, 0)
			}
		}
	}

	return
}
