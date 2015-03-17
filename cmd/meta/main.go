package main

import (
	"github.com/rakyll/command"
	"log"
	"os"
)

var (
	INFO  *log.Logger
	WARN  *log.Logger
	ERROR *log.Logger
)

func main() {
	// Register loggers.
	INFO = log.New(os.Stdout, "INFO: ", log.Ldate|log.Ltime|log.Lshortfile)
	WARN = log.New(os.Stdout, "WARN: ", log.Ldate|log.Ltime|log.Lshortfile)
	ERROR = log.New(os.Stderr, "ERROR:", log.Ldate|log.Ltime|log.Lshortfile)
	registerLogger()
	// Register commands.
	args := []string{}
	command.On("init", "generate strain information", &cmdInit{}, args)
	command.On("ortho_mcl", "find orthologs using OrthoMCL", &cmdOrthoMCL{}, args)
	command.On("ortho_aln", "align orthologs using MUSCLE", &cmdOrthoAln{}, args)
	command.On("cov_reads", "calculate correlation of subsitutions in reads", &cmdCovReads{}, args)
	command.On("cov_genomes", "calculate correlation of subsitutions in genomes", &cmdCovGenomes{}, args)
	command.On("bowtie2_index", "build bowtie2 index", &cmdIndex{}, []string{})
	command.On("bowtie2_align", "align reads using bowtie2", &cmdAlignReads{}, args)
	command.On("estimate", "estimate r/m", &cmdEstimate{}, args)

	// Parse and run commands.
	command.ParseAndRun()
}
