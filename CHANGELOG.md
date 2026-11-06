# Changelog

## 2.3.0

### Breaking changes

- Exit codes `200`, `201`, `202`, `203`, `204` and `222` are gone. Errors exit with `1` and the reason
  is logged at `error` level (#56).
- `eukcc folder --suffix` is removed. Bins are the files ending in `.fa`, `.fna` or `.fasta`, case
  insensitive and optionally gzip compressed (#55).
- `eukcc folder` stops with exit code `1` if any bin is not a valid FASTA file, empty files
  included. Every invalid file is logged, so all of them can be fixed before the next run.
- In `eukcc folder`, bins that could not be estimated (zero predicted proteins or no marker gene
  set) are listed in `eukcc.tsv` with `NA` and in `missing_marker_genes.txt`. Before, bins with
  zero proteins were in no output file.
- `eukcc define_set` exits with `1` when no marker gene set can be defined, as there is no result
  to write. The error lists the genomes and suggests a lower `--set_size` or `--marker_prevalence`.
- `eukcc single --gmes` is removed. It was legacy GeneMark-ES logic from EukCC 1 that no longer
  worked. Proteins are predicted with MetaEuk only.
- `eukcc.tsv` and `bad_quality.tsv` have a new column, `ncbi_lng_names`, after `ncbi_lng` (#51).
- `bad_quality.csv` and `merged_bins.csv` are renamed to `bad_quality.tsv` and `merged_bins.tsv`.
  They were already tab separated.
- Python 3.10 or 3.11 is required.
- `--debug` and `--quiet` are options of each subcommand and go after it
  (`eukcc single --debug ...`) (#40).
- `eukcc` uses `click`. `eukcc single --taxids` and `--genomes` take one value each and
  are repeated for several values (`--genomes a.fa --genomes b.fa`). All other options and
  defaults are the same.

### Added

- `ncbi_lng_names` column: the NCBI names of the taxids in `ncbi_lng`, in the same order,
  separated by `;` (#51).
- Gzip compressed (`.gz`) inputs, in `eukcc single` and `eukcc folder` (#41).
- Documentation of the database `eukcc2_db_ver_2026-11`.
- `eukcc download_db` downloads and extracts the database into `--out`. `--url` points to the
  current database on the EBI FTP by default and can be changed for other versions.

### Changed

- The sub-database (`fungi`, `protozoa`, `plant`) is chosen by a vote of the closest reference
  leaves of each placement, and needs at least 60% of the votes. It no longer depends on taxids
  in the code, so opisthokont protists can go to `protozoa`.
- The version is set only in `pyproject.toml` and read from the installed package.
  `eukcc/version.py` and `setup.py` are removed; the package is built with `python -m build`.
  The `binlinks.py` and `filter_euk_bins.py` commands, installed by `setup.py` until now, are
  kept as entry points, next to `binlinks` and `filter_euk_bins`.

### Fixed

- `eukcc ncbi_update` no longer fails looking for an output folder it does not have.
- The `filter_euk_bins` command works: the script had no `main` function for its entry point.

## [2.2.0](https://github.com/EBI-Metagenomics/EukCC/releases/tag/v.2.2.0) (2026-08-14)

### Fixed

- `tax_LCA` no longer fails with a `KeyError` when some placed proteins are not part of the marker
  gene set: they were taken as reference leaves, which have no taxid. Added a test (#54).

### Changed

- Quickstart documentation refers to the database version 1.2 and has typo fixes (#49).

## [2.1.3](https://github.com/EBI-Metagenomics/EukCC/releases/tag/v.2.1.3) (2025-01-31)

### Added

- `eukcc folder` writes `bad_quality.csv` and `no_placement_marker_genes.txt` (#46).
- Release to PyPI through a GitHub action (#47).

### Changed

- EukCC exits when MetaEuk predicts zero proteins (#46).
- Link to the newer database, logging, typo fixes and README updates, including the external
  tools needed by the PyPI and Bioconda installations (#46, #47).

## [2.0](https://github.com/EBI-Metagenomics/EukCC/releases/tag/v.2.0) (2021-06-09)

Major revision of EukCC. It breaks command line compatibility.

### Changed

- pplacer is replaced by EPA-ng.
- GeneMark-ES is replaced by MetaEuk.

## [0.3](https://github.com/EBI-Metagenomics/EukCC/releases/tag/v.0.3) (2021-02-01)

### Added

- Maximal silent contamination is reported when using GeneMark-ES.
- Users can provide a node to estimate completeness
  ([Finn-Lab/EukCC#15](https://github.com/Finn-Lab/EukCC/issues/15)).

### Changed

- The eukaryotic filter script is faster, easier to use and allows multithreading.

## [0.1.5.1](https://github.com/EBI-Metagenomics/EukCC/releases/tag/v.0.1.5.1) (2020-06-09)

### Changed

- Version bump for the Zenodo release of the paper, with the preprint information in the README.

## [0.1.5](https://github.com/EBI-Metagenomics/EukCC/releases/tag/v.0.1.5) (2020-06-08)

### Changed

- `eukcc.tsv` is created even if no estimation was done.
- Exit codes are more consistent.
- ete3 plotting routines are only imported when needed, so EukCC no longer depends on their
  installation.
