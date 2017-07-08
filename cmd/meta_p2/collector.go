package main

// CorrResult contains a correlation result.
type CorrResult struct {
	Lag      int
	Value    float64
	Count    int64
	Type     string
	Variance float64
}

// CorrResults is a list of CorrResult.
type CorrResults struct {
	GeneID  string
	Results []CorrResult
	ReadNum int
	GeneLen int
}

// Collector collect correlation results.
type Collector struct {
	m map[string][]*MeanVar
}

// NewCollector return a new Collector.
func NewCollector() *Collector {
	c := Collector{}
	c.m = make(map[string][]*MeanVar)
	return &c
}

// Add add an array of CorrResult.
func (c *Collector) Add(results CorrResults) {
	for _, res := range results.Results {
		for len(c.m[res.Type]) <= res.Lag {
			c.m[res.Type] = append(c.m[res.Type], NewMeanVar())
		}
		if res.Count > 0 {
			c.m[res.Type][res.Lag].Add(res.Value / float64(res.Count))
		}
	}
}

// Means return means of a particular type.
func (c *Collector) Means(corrType string) (values []float64) {
	for _, mv := range c.MeanVars(corrType) {
		values = append(values, mv.Mean())
	}
	return
}

// Vars return variances of a particular type.
func (c *Collector) Vars(corrType string) (values []float64) {
	for _, mv := range c.MeanVars(corrType) {
		values = append(values, mv.Variance())
	}
	return
}

// Ns return variances of a particular type.
func (c *Collector) Ns(corrType string) (nums []int) {
	for _, mv := range c.MeanVars(corrType) {
		nums = append(nums, mv.N)
	}
	return
}

// MeanVars return a list of meanvar.MeanVar.
func (c *Collector) MeanVars(corrType string) (values []*MeanVar) {
	return c.m[corrType]
}

// CorrTypes return all corr types.
func (c *Collector) CorrTypes() (corrTypes []string) {
	corrTypes = append(corrTypes, "P2")
	for key := range c.m {
		if key != "P2" {
			corrTypes = append(corrTypes, key)
		}
	}
	return
}

// Results get results
func (c *Collector) Results() (results []CorrResult) {
	corrTypes := c.CorrTypes()
	ks := 0.0
	for _, ctype := range corrTypes {
		means := c.Means(ctype)
		vars := c.Vars(ctype)
		ns := c.Ns(ctype)
		for i := 0; i < len(means); i++ {
			if ns[i] > 0 {
				res := CorrResult{}
				res.Lag = i * 3
				res.Count = int64(ns[i])
				res.Type = ctype
				res.Value = means[i]
				res.Variance = vars[i]
				if ctype == "P2" && i == 0 {
					res.Type = "Ks"
					ks = res.Value
				} else {
					if ks != 0 {
						res.Value /= ks
						res.Variance /= (ks * ks)
					}
				}
				results = append(results, res)
			}
		}
	}

	return
}
