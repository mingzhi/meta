package meta

// Strain information.
type Strain struct {
	Name        string   // taxonomy name.
	TaxId       string   // taxonomy Id.
	ProjectId   string   // project Id.
	Genomes     []Genome // genomes.
	Path        string   // path in referenc genome diretory.
	GeneticCode string   // genetic code id.
	Species     string   // species name
}

type Genome struct {
	Accession  string
	Replicon   string
	Length     int
	Seq        []byte
	PosProfile []byte
}
