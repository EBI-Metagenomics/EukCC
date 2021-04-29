# EukCC beta

![Coverage.py coverage](badges/coverage.svg)

EukCC is a completeness and contamination estimator for metagenomic assembled
and isolate genomes.

This is version 2, which should provide a better experience than
version 1. We aim at creating a stable package with long term support.

We are actuively working on this sofwtare and its not yet ready for release. Thats why the database is not yet 
availiable. EukCC2 will **not** work with EukCC1 databases.

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


