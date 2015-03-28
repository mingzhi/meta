package main

import (
	"encoding/json"
	"github.com/mingzhi/meta/ortho"
	"github.com/mingzhi/ncbiftp/seqrecord"
	"os"
	"path/filepath"
	"strings"
)

// Command to perform OrthoMCL similar functions.
type cmdOrthoMCL struct {
	cmdConfig // embed cmdConfig.
}

// Run command.
func (cmd *cmdOrthoMCL) Run(args []string) {
	// Parse config and settings.
	cmd.ParseConfig()
	cmd.LoadSpeciesMap()
	MakeDir(filepath.Join(*cmd.workspace, cmd.orthoOutBase))

	for prefix, strains := range cmd.speciesMap {
		if len(strains) >= 3 {
			INFO.Printf("%s\n", prefix)
			// OrthoMCL
			clusters := ortho.OrthoMCl(strains, cmd.refBase)

			// Write clusters into a file.
			cmd.writeClusters(prefix, clusters)

			// Find ortholog sequences.
			groups := ortho.FindOrthologs(strains, cmd.refBase, clusters)
			cmd.writeOrthologs(prefix, groups)
		} else {
			WARN.Printf("%s only has %d strains!\n", prefix, len(strains))
		}

	}

}

func (cmd *cmdOrthoMCL) writeClusters(prefix string, clusters [][]string) {
	fileName := prefix + ".mcl"
	filePath := filepath.Join(*cmd.workspace, cmd.orthoOutBase, fileName)
	f, err := os.Create(filePath)
	if err != nil {
		ERROR.Fatalln(err)
	}
	defer f.Close()

	for _, clst := range clusters {
		f.WriteString(strings.Join(clst, "\t") + "\n")
	}
}

func (cmd *cmdOrthoMCL) writeOrthologs(prefix string, groups []seqrecord.SeqRecords) {
	fileName := prefix + "_orthologs.json"
	filePath := filepath.Join(*cmd.workspace, cmd.orthoOutBase, fileName)
	w, err := os.Create(filePath)
	if err != nil {
		ERROR.Fatalln(err)
	}
	defer w.Close()

	encoder := json.NewEncoder(w)
	err = encoder.Encode(groups)
	if err != nil {
		ERROR.Fatalln(err)
	}
}
