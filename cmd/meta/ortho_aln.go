package main

import (
	"encoding/json"
	"github.com/mingzhi/meta/align/multi"
	"github.com/mingzhi/ncbiftp/seqrecord"
	"os"
	"path/filepath"
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

	for prefix, _ := range cmd.speciesMap {
		// Read ortholog groups.
		groups := cmd.ReadOrhtologs(prefix)

		if len(groups) > 0 {
			// Create a job for each sequence records.
			jobs := make(chan seqrecord.SeqRecords)
			go func() {
				defer close(jobs)
				for _, cluster := range groups {
					if len(cluster) >= 3 {
						jobs <- cluster
					}
				}
			}()

			// Create workers to do jobs.
			// done is signal channel.
			done := make(chan bool)
			// results is a channel for aligned sequence records.
			results := make(chan seqrecord.SeqRecords)
			for i := 0; i < *cmd.ncpu; i++ {
				go func() {
					for cluster := range jobs {
						aln := multi.Align(cluster, multi.Muscle)
						results <- aln
					}
					done <- true
				}()
			}

			// Waiting and checking done signal.
			go func() {
				defer close(results)
				for i := 0; i < *cmd.ncpu; i++ {
					<-done
				}
			}()

			// Collected aligned sequence records.
			alns := []seqrecord.SeqRecords{}
			for aln := range results {
				alns = append(alns, aln)
			}

			// Save aligned results into a json file.
			cmd.SaveAlignments(prefix, alns)
		} else {
			WARN.Printf("%s has zero orthologous cluster\n", prefix)
		}

	}

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

func (cmd *cmdOrthoAln) SaveAlignments(prefix string, alns []seqrecord.SeqRecords) {
	fileName := prefix + "_orthologs_aligned.json"
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
