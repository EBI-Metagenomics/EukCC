# EukCC 

![Coverage.py coverage](badges/coverage.svg)

EukCC is a completeness and contamination estimator for metagenomic assembled
microbial eukaryotic genomes.

## Documentation
Head over to https://eukcc.readthedocs.io/ to check out the documentation.

## Run

### Download EukCC2 database from [FTP](http://ftp.ebi.ac.uk/pub/databases/metagenomics/eukcc/)
```
# create a folder were to keep the database
mkdir eukccdb
cd eukccdb
wget http://ftp.ebi.ac.uk/pub/databases/metagenomics/eukcc/eukcc2_db_ver_2026-11.tar.gz
tar -xzvf eukcc2_db_ver_2026-11.tar.gz
export EUKCC2_DB=$(realpath eukcc2_db_ver_2026-11)
```

### Quickstart using container

Get EukCC quickly by fetching the BioContainers image:
https://quay.io/repository/biocontainers/eukcc

These images have no `latest` tag, so pick a version tag from the repository:

```
docker pull quay.io/biocontainers/eukcc:<tag>
singularity pull docker://quay.io/biocontainers/eukcc:<tag>
```

### Bioconda / pip

Alternatively you can install EukCC using conda or pip.

With pip, you also need to install these mandatory requirements:
- metaeuk=4.a0f584d
- pplacer
- epa-ng=0.3.8
- hmmer=3.3
- minimap2
- bwa
- samtools


## Outputs explanation

- `eukcc.log` - log of execution

### eukcc single

- `eukcc.tsv` - table with estimated completeness, contamination and taxonomy lineage (taxids and names)

### eukcc folder

- `eukcc.tsv` - table with estimated completeness, contamination and taxonomy lineage (taxids and names) for **good** quality bins
- `merged_bins.csv` - table of merged refined bins (tab separated)
- `bad_quality.csv` - tab separated table with estimated completeness, contamination and taxonomy lineage (taxids and names) for **bad** quality bins (chosen marker gene set is supported by less than half of the alignments)
- `missing_marker_genes.txt` - line separated list of bins with not defined set of marker genes
- `merged_bins` - folder with merged bins sequences 
- `refine_workdir` - working directory with intermediate steps results


## Don't use EukCC on already pubished data
Or at least not without thinking about it:

You should not use EukCC on already published genomes, if they were used during training of the marker
gene sets. If you want to make sure, the accessions used are listed in the `taxinfo` file of each
sub-database of the EukCC database (`files.taxinfo` in `db_<clade>/db_info.json`).

## Cite

If you use EukCC make sure to cite:

```
Saary, Paul, Alex L. Mitchell, and Robert D. Finn. 
"Estimating the quality of eukaryotic genomes recovered from metagenomic analysis with EukCC." 
Genome biology 21.1 (2020): 1-21.
```

EukCC also uses metaEUK, hmmer, pplacer, ete3 and epa-ng. 


## Changed compared to EukCC 1
**Note**: 
With version 2, EukCC should provide a better experience than
version 1. Version 2 is not compatible with previous versions, most commandline arguments changed.
So version 2 is **not** a drop in replacement.

- Users can set the prevalence threshold  for marker sets. In EukCC 1 
  this was fixed to 98% single copy prevalence.  Now users could change that to be more strict.
  We find that often 100% single copy prevalence can be found. 

## Issues and bugs

Please report any bugs and issues here on GitHub. Make sure to
include the debug log (run eukcc using `--debug` flag).

### used exit codes
- 0: Success.
- 1: Error, the reason is reported in the log
