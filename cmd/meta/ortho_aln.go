package main

import (
	"encoding/json"
	"github.com/mingzhi/biogo/seq"
	"github.com/mingzhi/meta/align/multi"
	"github.com/mingzhi/meta/genome"
	"github.com/mingzhi/meta/strain"
	"github.com/mingzhi/ncbiftp/seqrecord"
	"math"
	"os"
	"path/filepath"
	"sort"
	"strings"
)

// Command to align orthologs.
type cmdOrthoAln struct {
	cmdConfig // embed cmdConfig.
}

// Run command.
func (cmd *cmdOrthoAln) Run(args []string) {
	// Parse config and settings.
	cmd.ParseConfig()
	cmd.LoadSpeciesMap()
	MakeDir(filepath.Join(*cmd.workspace, cmd.orthoOutBase))

	for prefix, strains := range cmd.speciesMap {
		// Read ortholog protein clusters.
		rawClusters := cmd.ReadOrhtologs(prefix)
		// filter outliers based on their lengths.
		clusters := []seqrecord.SeqRecords{}
		for i := 0; i < len(rawClusters); i++ {
			records := rawClusters[i]
			if len(records) >= 3 {
				cls := filter(rawClusters[i])
				if len(cls) == len(records) {
					clusters = append(clusters, cls)
				}
			}

		}

		if len(clusters) > 0 {
			// align coding regions (protein clusters).
			alns := align(clusters, multi.AlignProt, *cmd.ncpu)
			cmd.SaveAlignments(prefix, alns)

			// expand gene to include its adjacent non-coding regions.
			appendix := "expanded"
			m := getGenomeMap(strains, cmd.refBase)
			expandedClusters := []seqrecord.SeqRecords{}
			for _, records := range clusters {
				expandedRecords := expand(records, m)
				if len(expandedRecords) > 0 {
					expandedClusters = append(expandedClusters, filter(expandedRecords))
				}
			}
			expandedAlns := align(expandedClusters, multi.AlignNucl, *cmd.ncpu)
			cmd.SaveAlignments(prefix, expandedAlns, appendix)
		} else {
			WARN.Printf("%s has zero orthologous cluster\n", prefix)
		}
	}
}

type multiAlignFunc func(seqRecords []seqrecord.SeqRecord, alignFunc multi.AlignFunc, options ...string) []seqrecord.SeqRecord

func align(clusters []seqrecord.SeqRecords, alignFunc multiAlignFunc, ncpu int) (alns []seqrecord.SeqRecords) {
	// Create a job for each sequence records.
	jobs := make(chan seqrecord.SeqRecords)
	go func() {
		defer close(jobs)
		for _, cluster := range clusters {
			if len(cluster) >= 3 {
				jobs <- cluster
			}
		}
	}()

	numWorker := ncpu

	// Create workers to do jobs.
	// done is signal channel.
	done := make(chan bool)
	// results is a channel for aligned sequence records.
	results := make(chan seqrecord.SeqRecords)
	for i := 0; i < numWorker; i++ {
		go func() {
			for cluster := range jobs {
				aln := alignFunc(cluster, multi.Muscle)
				results <- aln
			}
			done <- true
		}()
	}

	// Waiting and checking done signal.
	go func() {
		defer close(results)
		for i := 0; i < numWorker; i++ {
			<-done
		}
	}()

	// Collected aligned sequence records.
	for aln := range results {
		alns = append(alns, aln)
	}

	return
}

// return a map[string]genome.Genome
func getGenomeMap(strains []strain.Strain, refBase string) (genomeMap map[string]genome.Genome) {
	genomeMap = make(map[string]genome.Genome)
	genomes := loadGenomes(strains, refBase)
	for _, g := range genomes {
		genomeMap[g.RefAcc()] = g
	}
	return
}

// for each genome in a strain, load its DNA sequence and position profile.
func loadGenomes(strains []strain.Strain, refBase string) []genome.Genome {
	genomes := []genome.Genome{}
	for _, s := range strains {
		base := filepath.Join(refBase, s.Path)
		for _, g := range s.Genomes {
			genome.LoadFna(&g, base)
			genome.LoadProfile(&g, base)
			genomes = append(genomes, g)
		}
	}

	return genomes
}

// increasing index by specific direction -- plus or minor
func updateIndex(i int, plus bool) int {
	if plus {
		return i + 1
	} else {
		return i - 1
	}
}

const (
	AlphabetDNA = "ATGCatgc"
)

func isValidNucl(r byte) bool {
	for i := 0; i < len(AlphabetDNA); i++ {
		if AlphabetDNA[i] == r {
			return true
		}
	}
	return false
}

func expand(originals seqrecord.SeqRecords, m map[string]genome.Genome) seqrecord.SeqRecords {
	updates := seqrecord.SeqRecords{}
	for _, r := range originals {
		g, found := m[genome.FindRefAcc(r.Genome)]
		if found {
			s := expandOne(r, g)
			updates = append(updates, s)
		} else {
			WARN.Printf("Can not find genome for %s, with genome accession %s\n", r.Id, r.Genome)
		}
	}
	return updates
}

func expandOne(r seqrecord.SeqRecord, g genome.Genome) (s seqrecord.SeqRecord) {
	directions := []bool{false, true}
	ends := []int{}
	for i := 0; i < len(directions); i++ {
		plus := directions[i]

		start := r.Loc.From
		if plus {
			start = r.Loc.To
		}
		start = start - 1
		for {
			start = updateIndex(start, plus)
			index := (len(g.Seq) + start) % len(g.Seq)
			pos := g.PosProfile[index]
			if pos == 0 {
				break
			}
		}

		for {
			index := (len(g.Seq) + start) % len(g.Seq)
			pos := g.PosProfile[index]
			nuc := g.Seq[index]
			if pos != 0 || !isValidNucl(nuc) {
				break
			}
			start = updateIndex(start, plus)
		}

		ends = append(ends, start)
	}

	sort.Ints(ends)
	length := len(g.Seq)
	from, to := (ends[0]+1+length)%length, (ends[1]+length)%length
	nucl := []byte{}
	if from > to {
		nucl = g.Seq[from:]
		nucl = g.Seq[:to]
	} else {
		nucl = g.Seq[from:to]
	}

	s.Code = r.Code
	s.Genome = r.Genome
	s.Id = r.Id
	s.Name = r.Name
	s.Loc.From = from + 1
	s.Loc.To = to
	s.Loc.Strand = r.Loc.Strand
	if s.Loc.Strand == "-" {
		nucl = seq.Complement(seq.Reverse(nucl))
	}
	s.Nucl = nucl

	return s
}

// filter outliers based on the length to median.
func filter(records seqrecord.SeqRecords) seqrecord.SeqRecords {
	lengths := []int{}
	for i := 0; i < len(records); i++ {
		lengths = append(lengths, len(records[i].Nucl))
	}
	sort.Ints(lengths)

	l := len(lengths)
	median := (lengths[l/2] + lengths[(l-1)/2]) / 2
	news := seqrecord.SeqRecords{}
	for i := 0; i < len(records); i++ {
		ratio := float64(len(records[i].Nucl)-median) / float64(median)
		if math.Abs(ratio) <= 0.1 {
			news = append(news, records[i])
		}
	}
	return news
}

func (cmd *cmdOrthoAln) ReadOrhtologs(prefix string) (groups []seqrecord.SeqRecords) {
	fileName := prefix + "_orthologs.json"
	filePath := filepath.Join(*cmd.workspace, cmd.orthoOutBase,
		fileName)
	r, err := os.Open(filePath)
	if err != nil {
		WARN.Println(err)
		return
	}
	defer r.Close()

	decoder := json.NewDecoder(r)
	err = decoder.Decode(&groups)
	if err != nil {
		ERROR.Fatalln(err)
	}

	return
}

func (cmd *cmdOrthoAln) SaveAlignments(prefix string, alns []seqrecord.SeqRecords, appendix ...string) {
	prefixTerms := []string{prefix}
	if len(appendix) > 0 {
		prefixTerms = append(prefixTerms, appendix...)
	}

	fileName := strings.Join(prefixTerms, "_") + "_orthologs_aligned.json"
	filePath := filepath.Join(*cmd.workspace, cmd.orthoOutBase,
		fileName)
	w, err := os.Create(filePath)
	if err != nil {
		ERROR.Fatalln(err)
	}
	defer w.Close()

	encoder := json.NewEncoder(w)
	err = encoder.Encode(alns)
	if err != nil {
		ERROR.Fatalln(err)
	}
}
