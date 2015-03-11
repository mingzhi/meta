package main

import (
	"encoding/json"
	"github.com/mingzhi/meta"
	"github.com/mingzhi/ncbiutils"
	"log"
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
	MakeDir(cmd.orthoOutBase)

	// Read strain information.
	strainFilePath := filepath.Join(*cmd.workspace, cmd.strainFileName)
	strains := meta.ReadStrains(strainFilePath)

	// OrthoMCL
	clusters := meta.OrthoMCl(strains, cmd.refBase)

	// Write clusters into a file.
	cmd.writeClusters(clusters)

	// Find ortholog sequences.
	groups := meta.FindOrthologs(strains, cmd.refBase, clusters)
	cmd.writeOrthologs(groups)
}

func (cmd *cmdOrthoMCL) writeClusters(clusters [][]string) {
	fileName := cmd.prefix + ".mcl"
	filePath := filepath.Join(*cmd.workspace, cmd.orthoOutBase, fileName)
	f, err := os.Create(filePath)
	if err != nil {
		log.Panic(err)
	}
	defer f.Close()

	for _, clst := range clusters {
		f.WriteString(strings.Join(clst, "\t") + "\n")
	}
}

func (cmd *cmdOrthoMCL) writeOrthologs(groups []ncbiutils.SeqRecords) {
	fileName := cmd.prefix + "_orthologs.json"
	filePath := filepath.Join(*cmd.workspace, cmd.orthoOutBase, fileName)
	w, err := os.Create(filePath)
	if err != nil {
		log.Panic(err)
	}
	defer w.Close()

	encoder := json.NewEncoder(w)
	err = encoder.Encode(groups)
	if err != nil {
		log.Panic(err)
	}
}
