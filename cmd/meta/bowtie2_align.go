package main

import (
	"bytes"
	"fmt"
	"github.com/mingzhi/meta"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"strings"
)

var (
	bowtiedSamAppendix string = ".bowtie2_aligned.sam"
	bowtiedLogAppendix string = ".bowtie2_aligned.log"
)

// Command for mapping reads to reference genomes.
type cmdAlignReads struct {
	cmdConfig // embedded cmdConfig.
}

func (cmd *cmdAlignReads) Run(args []string) {
	// Parse configure and settings.
	cmd.ParseConfig()
	MakeDir(filepath.Join(*cmd.workspace, cmd.samOutBase))

	// Map reads to each strain.
	jobs := make(chan meta.Strain)
	go func() {
		defer close(jobs)
		for _, strains := range cmd.speciesMap {
			for _, s := range strains {
				jobs <- s
			}
		}
	}()

	// Create cmd.ncpu workers for aligning.
	// send done signal when the job is done.
	done := make(chan bool)
	for i := 0; i < *cmd.ncpu; i++ {
		go func() {
			for strain := range jobs {
				cmd.align(strain)
			}
			done <- true
		}()
	}

	// Waiting for workers.
	for i := 0; i < *cmd.ncpu; i++ {
		<-done
	}
}

// bowie2 options:
//  -q:fastq, -f:fasta
//  -1 <paired-end-1> -2 <paired-end-2> or -U <single-end>
//  -p %s : Launch NTHREADS parallel search threads (default: 1)
//  -t : Print the wall-clock time
//  --no-unal : Suppress SAM records for reads that failed to align.
//  -I %s : The minimum fragment length for valid paired-end alignments (default: 0)
//  -X %s : The maximum fragment length for valid paired-end alignments (default: 500)
//  --no-discordant : A discordant alignment is an alignment where both mates align
//                    uniquely, but that does not satisfy the paired-end constraints.
//                    --no-discordant only report concordant alignments.
//  --no-mixed : By default, when bowtie2 cannot find a concordant or discordant
//               alignment for a pair, it then tries to find alignments for the individual
//               mates. This option disables that behavior.
//  --ignore-quals : Ignore qualities of mismatches
//  -mp 6,6 : Sets the maximum (MX) and minimum (MN) mismatch penalties to 6
//  -np 6 : Sets penalty for positions where the read, reference, or both, contain an
//          ambiguous character such as N. Default: 1.
//  --score-min <func> : Sets a function governing the minimum alignment score needed for
//                       an alignment to be considered "valid"
//                       L,0,-0.6 sets the minimum-score function f to f(x) = 0 + -0.6 * x,
//                       where x is the read length. For example, if we set L,-18.0,0.0,
//                       then it allows "3 (18/6)" mismatches without considering length of read.
//  --gbar <int> : Disallow gaps within <int> positions of the beginning or end of the read.
//                 Default: 4. We set as 1000 for not allowing gaps.
//
func (cmd *cmdAlignReads) align(strain meta.Strain) {
	// Basic control options.
	options := []string{
		"-t",
		"--no-unal",
		"--no-discordant",
		"--no-mixed",
		"--ignore-quals",
	}

	// for fasta or fastq
	// here default to be fastq.
	options = append(options, "-q")
	// bowtie threads number.
	threadsNum := strconv.Itoa(cmd.bowtieThreadsNum)
	options = append(options, []string{"-p", threadsNum}...)

	genomeIndexBase := filepath.Join(cmd.refBase, strain.Path, strain.Path)
	outFilePrefix := filepath.Join(*cmd.workspace, cmd.samOutBase, strain.Path)
	samOutFilePath := outFilePrefix + bowtiedSamAppendix
	options = append(options, []string{"-x", genomeIndexBase}...)
	options = append(options, []string{"-1", cmd.pairedEndReadFile1}...)
	options = append(options, []string{"-2", cmd.pairedEndReadFile2}...)
	options = append(options, []string{"-S", samOutFilePath}...)

	// calculate min_score based on number of mismatches.
	if cmd.maximumMismatchCount > 0 {
		min_score := cmd.maximumMismatchCount * (-6)
		options = append(options, []string{"--mp", "6,6"}...)
		options = append(options, []string{"--np", "6"}...)
		options = append(options, []string{"--score-min", fmt.Sprintf("L,%d,0.0", min_score)}...)
		options = append(options, []string{"--gbar", "1000"}...)
	}

	command := exec.Command("bowtie2", options...)
	stderr := new(bytes.Buffer)
	command.Stderr = stderr
	logFilePath := outFilePrefix + bowtiedLogAppendix
	logFile, err := os.Create(logFilePath)
	if err != nil {
		ERROR.Printf("Cannot create %s: %v\n", logFilePath, err)
	}
	defer logFile.Close()
	command.Stdout = logFile
	err = command.Run()
	if err != nil {
		ERROR.Println(strings.Join(options, " "))
		ERROR.Println(err)
		ERROR.Fatalln(string(stderr.Bytes()))
	}
}
