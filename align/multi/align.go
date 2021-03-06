package multi

import (
	"bytes"
	"github.com/mingzhi/biogo/seq"
	"github.com/mingzhi/ncbiftp/seqrecord"
	"io"
	"os/exec"
	"strings"
)

type AlignFunc func(stdin io.Reader, stdout, stderr io.Writer, options ...string) error

// Multiple sequence alignment of protein sequences
// and back translate them to nucleotide sequences
func AlignProt(seqRecords []seqrecord.SeqRecord, alignFunc AlignFunc, options ...string) []seqrecord.SeqRecord {
	// prepare protein sequences in fasta format
	stdin := new(bytes.Buffer)
	srMap := make(map[string]seqrecord.SeqRecord)
	for _, sr := range seqRecords {
		stdin.WriteString(">" + sr.Id + "|" + sr.Genome + "\n")
		stdin.Write(sr.Prot)
		stdin.WriteString("\n")
		srMap[sr.Id+"|"+sr.Genome] = sr
	}
	stdout := new(bytes.Buffer)
	stderr := new(bytes.Buffer)
	alignFunc(stdin, stdout, stderr, options...)
	fr := seq.NewFastaReader(stdout)
	alns, err := fr.ReadAll()
	if err != nil {
		panic(string(stderr.Bytes()))
	}

	alnSeqRecords := []seqrecord.SeqRecord{}
	for _, aaSeq := range alns {
		if aaSeq != nil {
			if _, found := srMap[aaSeq.Id]; found {
				na := BackTranslate(aaSeq.Seq, srMap[aaSeq.Id].Nucl)
				sr := seqrecord.SeqRecord{
					Id:     strings.Split(aaSeq.Id, "|")[0],
					Prot:   aaSeq.Seq,
					Nucl:   na,
					Genome: srMap[aaSeq.Id].Genome,
					Name:   srMap[aaSeq.Id].Name,
					Loc:    srMap[aaSeq.Id].Loc,
					Code:   srMap[aaSeq.Id].Code,
				}
				alnSeqRecords = append(alnSeqRecords, sr)
			}
		}
	}

	return alnSeqRecords

}

// Multiple sequence alignment of protein sequences
// and back translate them to nucleotide sequences
func AlignNucl(seqRecords []seqrecord.SeqRecord, alignFunc AlignFunc, options ...string) []seqrecord.SeqRecord {
	// prepare protein sequences in fasta format
	stdin := new(bytes.Buffer)
	srMap := make(map[string]seqrecord.SeqRecord)
	for _, sr := range seqRecords {
		stdin.WriteString(">" + sr.Id + "|" + sr.Genome + "\n")
		stdin.Write(sr.Nucl)
		stdin.WriteString("\n")
		srMap[sr.Id+"|"+sr.Genome] = sr
	}
	stdout := new(bytes.Buffer)
	stderr := new(bytes.Buffer)
	alignFunc(stdin, stdout, stderr, options...)
	fr := seq.NewFastaReader(stdout)
	alns, err := fr.ReadAll()
	if err != nil {
		panic(string(stderr.Bytes()))
	}

	alnSeqRecords := []seqrecord.SeqRecord{}
	for _, aaSeq := range alns {
		if aaSeq != nil {
			if _, found := srMap[aaSeq.Id]; found {
				sr := seqrecord.SeqRecord{
					Id:     strings.Split(aaSeq.Id, "|")[0],
					Nucl:   aaSeq.Seq,
					Prot:   srMap[aaSeq.Id].Prot,
					Genome: srMap[aaSeq.Id].Genome,
					Name:   srMap[aaSeq.Id].Name,
					Loc:    srMap[aaSeq.Id].Loc,
					Code:   srMap[aaSeq.Id].Code,
				}
				alnSeqRecords = append(alnSeqRecords, sr)
			}
		}
	}

	return alnSeqRecords

}

// do multiple sequence alignment using muscle
func Muscle(stdin io.Reader, stdout, stderr io.Writer, options ...string) (err error) {
	cmd := exec.Command("muscle", options...)
	cmd.Stdin = stdin
	cmd.Stdout = stdout
	cmd.Stderr = stderr
	err = cmd.Run()
	return
}

// back translate amino acid alignment to nucleotide sequences.
func BackTranslate(aa, na []byte) []byte {
	k := 0
	aln := []byte{}
	for i := 0; i < len(aa); i++ {
		if aa[i] == '-' {
			aln = append(aln, []byte{'-', '-', '-'}...)
		} else {
			aln = append(aln, na[k*3:(k+1)*3]...)
			k++
		}
	}
	return aln
}
