package main

import (
	"encoding/json"
	"flag"
	"github.com/mingzhi/meta"
	"github.com/mingzhi/ncbiutils"
	"github.com/rakyll/command"
	"log"
	"os"
	"path/filepath"
	"regexp"
	"runtime"
	"strings"
)

var (
	DefaultMaxProcs = runtime.NumCPU()
	INFO            *log.Logger
	WARN            *log.Logger
)

func main() {
	runtime.GOMAXPROCS(DefaultMaxProcs)
	// Register loggers.
	INFO = log.New(os.Stdout, "INFO: ", log.Ldate|log.Ltime|log.Lshortfile)
	WARN = log.New(os.Stdout, "WARN: ", log.Ldate|log.Ltime|log.Lshortfile)
	// Register commands.
	command.On("init", "initialize", &cmdInit{}, []string{})
	command.On("makeudb", "make usearch db", &makeUdbCmd{}, []string{})
	command.On("orthomcl", "OrthoMCL", &ortholMclCmd{}, []string{})
	command.On("orthoaln", "align orthologs", &alignOrthologCmd{}, []string{})
	command.On("readanchor", "read anchor", &readAnchorCmd{}, []string{})
	command.On("profile", "calculate position profile for a referenc genome",
		&cmdProfile{}, []string{})
	command.On("cov", "calculate cov for reads mapped to a reference genome",
		&cmdCov{}, []string{})
	// Parse and run commands.
	command.ParseAndRun()
}

func registerLogger() {
	meta.Info = log.New(os.Stdout, "INFO: ", log.Ldate|log.Ltime|log.Lshortfile)
	meta.Warn = log.New(os.Stdout, "WARN: ", log.Ldate|log.Ltime|log.Lshortfile)
}

type ortholMclCmd struct {
	workspace *string
	ref       *string
	prefix    *string
}

func (cmd *ortholMclCmd) Flags(fs *flag.FlagSet) *flag.FlagSet {
	cmd.workspace = fs.String("w", "", "workspace")
	cmd.ref = fs.String("r", "", "reference genome diretory")
	cmd.prefix = fs.String("p", "", "prefix")

	return fs
}

func (cmd *ortholMclCmd) Run(args []string) {
	registerLogger()
	speciesMapFileName := filepath.Join(*cmd.workspace, "species_map.json")
	speciesMap := meta.ReadSpeciesMap(speciesMapFileName)
	strains, found := speciesMap[*cmd.prefix]
	if !found {
		log.Fatalf("Can not find strain information for %s from species_map.json\n", *cmd.prefix)
	}
	clusters := meta.OrthoMCl(strains, *cmd.ref)

	// Write the clusters into a file.
	filePath := filepath.Join(*cmd.workspace, *cmd.prefix+".mcl")
	w, err := os.Create(filePath)
	if err != nil {
		log.Panic(err)
	}
	defer w.Close()
	for _, clt := range clusters {
		w.WriteString(strings.Join(clt, "\t") + "\n")
	}

	INFO.Printf("OrthoMCL clusters were saved to %s\n", filePath)

	orthologGroups := meta.FindOrthologs(strains, *cmd.ref, clusters)
	filePath2 := filepath.Join(*cmd.workspace, *cmd.prefix+"_orthologs.json")
	w2, err := os.Create(filePath2)
	if err != nil {
		log.Panic(err)
	}
	defer w2.Close()
	ec := json.NewEncoder(w2)
	err = ec.Encode(orthologGroups)
	if err != nil {
		log.Panic(err)
	}
}

type makeUdbCmd struct {
	workspace *string
	ref       *string
	prefix    *string
}

func (cmd *makeUdbCmd) Flags(fs *flag.FlagSet) *flag.FlagSet {
	cmd.workspace = fs.String("w", "", "workspace")
	cmd.ref = fs.String("r", "", "reference genome diretory")
	cmd.prefix = fs.String("p", "", "prefix")

	return fs
}

func (cmd *makeUdbCmd) Run(args []string) {
	registerLogger()
	speciesMapFileName := filepath.Join(*cmd.workspace, "species_map.json")
	speciesMap := meta.ReadSpeciesMap(speciesMapFileName)
	strains, found := speciesMap[*cmd.prefix]
	if !found {
		log.Fatalf("Can not find strain information for %s from species_map.json\n", *cmd.prefix)
	}

	for _, s := range strains {
		for _, g := range s.Genomes {
			acc := cleanAccession(g.Accession)
			fileName := filepath.Join(*cmd.ref, s.Path, acc+".faa")
			meta.UsearchMakeUDB(fileName)
		}
	}
}

func cleanAccession(name string) string {
	return strings.Split(name, ".")[0]
}

type alignOrthologCmd struct {
	workspace *string
	prefix    *string
}

func (cmd *alignOrthologCmd) Flags(fs *flag.FlagSet) *flag.FlagSet {
	cmd.workspace = fs.String("w", "", "workspace")
	cmd.prefix = fs.String("p", "", "prefix")
	return fs
}

func (cmd *alignOrthologCmd) Run(args []string) {
	// Read ortholog clusters.
	filePath := filepath.Join(*cmd.workspace, *cmd.prefix+"_orthologs.json")
	f, err := os.Open(filePath)
	if err != nil {
		log.Panic(err)
	}
	defer f.Close()

	dc := json.NewDecoder(f)
	orthologs := []ncbiutils.SeqRecords{}
	err = dc.Decode(&orthologs)
	if err != nil {
		panic(err)
	}

	ncpu := runtime.GOMAXPROCS(0)
	jobs := make(chan ncbiutils.SeqRecords)
	go func() {
		for _, cluster := range orthologs {
			if len(cluster) >= 3 {
				jobs <- cluster
			}
		}
		close(jobs)
	}()

	done := make(chan bool)
	results := make(chan ncbiutils.SeqRecords)
	for i := 0; i < ncpu; i++ {
		go func() {
			for cluster := range jobs {
				aln := meta.MultiAlign(cluster, meta.Muscle)
				results <- aln
			}
			done <- true
		}()
	}

	go func() {
		for i := 0; i < ncpu; i++ {
			<-done
		}
		close(results)
	}()

	alignments := []ncbiutils.SeqRecords{}
	for aln := range results {
		alignments = append(alignments, aln)
	}

	outPath := filepath.Join(*cmd.workspace, *cmd.prefix+"_orthologs_aln.json")
	w, err := os.Create(outPath)
	if err != nil {
		log.Panic(err)
	}
	defer w.Close()

	ec := json.NewEncoder(w)
	err = ec.Encode(alignments)
	if err != nil {
		log.Panic(err)
	}
}

type readAnchorCmd struct {
	workspace *string
	prefix    *string
	samout    *string
}

func (cmd *readAnchorCmd) Flags(fs *flag.FlagSet) *flag.FlagSet {
	cmd.workspace = fs.String("w", "", "workspace")
	cmd.prefix = fs.String("p", "", "prefix")
	cmd.samout = fs.String("s", "", "bam output folder")

	return fs
}

func (cmd *readAnchorCmd) Run(args []string) {
	registerLogger()
	// find strains.
	speciesMapFileName := filepath.Join(*cmd.workspace, "species_map.json")
	speciesMap := meta.ReadSpeciesMap(speciesMapFileName)
	strains, found := speciesMap[*cmd.prefix]
	if !found {
		WARN.Fatalf("Can't not find strains for %s\n", *cmd.prefix)
	}

	readMap := make(map[string]meta.SamRecords)
	for _, strain := range strains {
		fileName := strain.Path + ".align.bam"
		filePath := filepath.Join(*cmd.samout, strain.Path, fileName)
		header, records := meta.ReadBamFile(filePath)
		m := meta.SeparateSamRecords(header.Refs(), records)
		for genome, recs := range m {
			acc := findRefAcc(genome)
			readMap[acc] = recs
			INFO.Println(genome)
		}
	}

	alignments := meta.ReadAlignments(filepath.Join(*cmd.workspace, *cmd.prefix+"_orthologs_aln.json"))

	mappedReads := meta.ReadAnchor(readMap, alignments)
	n := 0
	for _, mr := range mappedReads {
		if len(mr) > 500 {
			n++
		}
	}
	INFO.Println(n)
}

func findRefAcc(name string) string {
	re := regexp.MustCompile("NC_\\d+")
	return re.FindString(name)
}
