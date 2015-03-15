package main

import (
	"encoding/json"
	"github.com/mingzhi/meta"
	"github.com/mingzhi/ncbiutils"
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
	MakeDir(filepath.Join(*cmd.workspace, cmd.orthoOutBase))

	for prefix, strains := range cmd.speciesMap {
		INFO.Printf("%s\n", prefix)
		// OrthoMCL
		clusters := meta.OrthoMCl(strains, cmd.refBase)

		// Write clusters into a file.
		cmd.writeClusters(prefix, clusters)

		// Find ortholog sequences.
		groups := meta.FindOrthologs(strains, cmd.refBase, clusters)
		cmd.writeOrthologs(prefix, groups)
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

func (cmd *cmdOrthoMCL) writeOrthologs(prefix string, groups []ncbiutils.SeqRecords) {
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
