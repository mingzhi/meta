package reads

import (
	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
	"io"
	"log"
	"os"
)

// Read SAM file and return its header and records.
// NOT explicitly sorted.
func ReadSamFile(fileName string) (header *sam.Header, records []*sam.Record) {
	// Open sam file.
	f, err := os.Open(fileName)
	if err != nil {
		log.Panic(err)
	}
	defer f.Close()

	// Create sam reader,
	// and read the reference genomes.
	reader, err := sam.NewReader(f)
	if err != nil {
		log.Panic(err)
	}
	header = reader.Header()

	for {
		r, err := reader.Read()
		if err != nil {
			if err != io.EOF {
				log.Panic(err)
			} else {
				break
			}
		}
		records = append(records, r)
	}
	return
}

// Read BAM file and return its header and records.
// NOT explicitly sorted.
func ReadBamFile(fileName string) (header *sam.Header, records []*sam.Record) {
	// Open bam file.
	f, err := os.Open(fileName)
	if err != nil {
		log.Panic(err)
	}
	defer f.Close()

	// Create bam reader,
	// and read the reference genomes.
	rd := 0 // ignore this now.
	reader, err := bam.NewReader(f, rd)
	if err != nil {
		log.Panic(err)
	}
	header = reader.Header()

	for {
		r, err := reader.Read()
		if err != nil {
			if err != io.EOF {
				log.Panic(err)
			}
			break
		}
		records = append(records, r)
	}

	return
}
