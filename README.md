# EukCC 

![Coverage.py coverage](badges/coverage.svg)

EukCC is a completeness and contamination estimator for metagenomic assembled
and isolate genomes.

This is version 2, which should provide a better experience than
version 1. We aim at creating a stable package with long term support.

## Documentation
Got over to eukcc.readthedocs.io/ to check out the documentation.

## Get the container

Get EukCC quickly by fetching the container.
```
docker pull ...
singularity pull docker://
```
Make sure to also fetch the database for version 2:

http://ftp.ebi.ac.uk/pub/databases/metagenomics/eukcc/

## Cite

If you use EukCC make sure to cite:

```
Saary, Paul, Alex L. Mitchell, and Robert D. Finn. "Estimating the quality of eukaryotic genomes recovered from metagenomic analysis with EukCC." Genome biology 21.1 (2020): 1-21.
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


