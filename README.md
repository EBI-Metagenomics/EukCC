# EukCC 

![Coverage.py coverage](badges/coverage.svg)

EukCC is a completeness and contamination estimator for metagenomic assembled
microbial eukaryotic genomes. 

With version 2, EukCC should provide a better experience than
version 1. We aim at creating a stable package with long term support.
Version 2 is currently under development and using it could lead to unexpected issues. We aim to 
release EukCC version 2 soon.

**Note**: Version 2 is not compatible with previous versions, most commandline arguments changed.
So version 2 is **not** a drop in replacement.

## Documentation
Head over to https://eukcc.readthedocs.io/ to check out the documentation.


## Quickstart

Here a super fast overview on how to get started using singularity. For more information
please see the documentation.

```
mkdir eukccdb
cd eukccdb
wget http://ftp.ebi.ac.uk/pub/databases/metagenomics/eukcc/eukcc2_db_ver_1.1.tar.gz
tar -xzvf eukcc2_db_ver_1.1.tar.gz
export EUKCC2_DB=$(realpath eukcc2_db_ver_1.1)
singularity pull docker://microbiomeinformatics/eukcc

singularity exec eukcc2_latest.sif eukcc single -h
```

## Get the container

Get EukCC quickly by fetching the container [Container repro will change soon]
```
docker pull microbiomeinformatics/eukcc
singularity pull docker://microbiomeinformatics/eukcc
```
If you installed version 2: make sure to also fetch the database for version 2:

http://ftp.ebi.ac.uk/pub/databases/metagenomics/eukcc/


## Don't use EukCC on already pubished data
Or at least not without thinking about it:

You should not use EukCC on already published genomes, if they have used during training of the marker 
gene sets. If you want to make sure, you can see all used accessions in the database file `db_base/backbone/base_taxinfo.csv`.

## Cite

If you use EukCC make sure to cite:

```
Saary, Paul, Alex L. Mitchell, and Robert D. Finn. 
"Estimating the quality of eukaryotic genomes recovered from metagenomic analysis with EukCC." 
Genome biology 21.1 (2020): 1-21.
```

EukCC also uses metaEUK, hmmer, pplacer, ete3 and epa-ng. 


## Changed compared to EukCC 1
- Users can set the prevalence threshold  for marker sets. In EukCC 1 
  this was fixed to 98% single copy prevalence.  Now users could change that to be more strict.
  We find that often 100% single copy prevalence can be found. 

## Issues and bugs

Please report any bugs and issues here on GitHub. Make sure to
include the debug log (run eukcc using `--debug` flag).

### used exit codes
- 200: File not found
- 201: No Marker gene set could be defined
- 202: No database provided
- 203: Corrupted file
- 222: Invalid settings


