package main

import (
	"bytes"
	"flag"
	"github.com/mingzhi/meta"
	"github.com/spf13/viper"
	"log"
	"os/exec"
	"path/filepath"
	"runtime"
	"strconv"
)

// Command for mapping reads to reference genomes.
type cmdAlignReads struct {
	workspace *string // workspace.
	config    *string // configure file name.

	refBase            string // reference genome folder.
	strainFileName     string // strain file name.
	pairedEndReadFile1 string // paired-end read file 1.
	pairedEndReadFile2 string // paired-end read file 2.
	ncpu               int    // number of CPUs for using.
	bowtieThreadsNum   int    // bowtie threads number.
	samOutBase         string // sam output folder.
}

// Register flags.
func (cmd *cmdAlignReads) Flags(fs *flag.FlagSet) *flag.FlagSet {
	cmd.workspace = fs.String("w", "", "workspace")
	cmd.config = fs.String("c", "config", "configure file name")
	return fs
}

// Initialize.
func (cmd *cmdAlignReads) Init() {
	// Register viper for configurations.
	viper.SetConfigName(*cmd.config)
	viper.AddConfigPath(*cmd.workspace)
	viper.ReadInConfig()

	// Read settings.
	cmd.refBase = viper.GetString("Reference_Genome_Directory")
	cmd.strainFileName = viper.GetString("Strain_File_Name")
	cmd.pairedEndReadFile1 = viper.GetString("Paired_End_Reads_1")
	cmd.pairedEndReadFile2 = viper.GetString("Paired_End_Reads_2")
	cmd.bowtieThreadsNum = viper.GetInt("Bowtie_Threads_Num")
	cmd.samOutBase = viper.GetString("Sam_Output_Diretory")
	if cmd.bowtieThreadsNum <= 0 {
		cmd.bowtieThreadsNum = 1
	}
	cmd.ncpu = runtime.GOMAXPROCS(0)
}

func (cmd *cmdAlignReads) Run(args []string) {
	cmd.Init()
	// Read strain information.
	strainFilePath := filepath.Join(*cmd.workspace, cmd.strainFileName)
	strains := meta.ReadStrains(strainFilePath)
	// Map reads to each strain.
	jobs := make(chan meta.Strain)
	go func() {
		for _, s := range strains {
			jobs <- s
		}
		close(jobs)
	}()

	// Create cmd.ncpu workers for aligning.
	// send done signal when the job is done.
	done := make(chan bool)
	for i := 0; i < cmd.ncpu; i++ {
		go func() {
			for strain := range jobs {
				cmd.align(strain)
			}
			done <- true
		}()
	}

	// Waiting for workers.
	for i := 0; i < cmd.ncpu; i++ {
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
	samOutFilePath := filepath.Join(cmd.samOutBase, strain.Path+".sam")
	options = append(options, []string{"-x", genomeIndexBase}...)
	options = append(options, []string{"-1", cmd.pairedEndReadFile1}...)
	options = append(options, []string{"-2", cmd.pairedEndReadFile2}...)
	options = append(options, []string{"-S", samOutFilePath}...)
	command := exec.Command("bowtie2", options...)
	stderr := new(bytes.Buffer)
	stdout := new(bytes.Buffer)
	command.Stderr = stderr
	command.Stdout = stdout
	err := command.Run()
	if err != nil {
		log.Println(string(stdout.Bytes()))
		log.Fatalln(string(stderr.Bytes()))
	}
}
