
sample_info.csv contains information about each genotyped sample
	Name: sample name used in this dataset
	Population: either "BSCB" or "BSSS"
	Type: whether the sample is from one of the populations, or is a founder / derived line of one of the populations
	Cycle: The cycle of selection from which the sample originates, or from which the line was derived.

SNP_info.csv contains information about each of the genotyped SNPs
	Name: The SNP name, this is used in the genotype_data.csv file
	ID: The SNP ID on the 56k platform
	Chromosome: The chromosome to which the SNP was mapped
	Position: The physical position used in this study
	Genetic_Position: The genetic position used in this study, in centimorgans


genotype_data.csv contains the actual genotype calls.
	The first row contains the sample names.
	Each following row contains the marker name as the first entry, and the genotypes of the samples in the following entries.
