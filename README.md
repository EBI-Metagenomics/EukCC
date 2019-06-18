# EukCC

Estimate completeness and contamination in Eukaryotes.

**this software is under development and unstable. Updates might break or 
enhance your experience**

## Install

### Dependencies

**GeneMark-ES**

If you want to predict genes install GeneMark-ES: http://exon.gatech.edu/GeneMark/license_download.cgi

For this you need some perl libraries installed:

```
cpan install YAML
cpan install Logger::Simple
cpan install Parallel::ForkManager.pm
```

**HMMER**

Just install HMMER from: http://hmmer.org/

I used version 3.2.1 but any later version should work fine.


**Pplacer**

Install pplacer. I installed it via conda:

```
conda install -c bioconda pplacer 
```

I used version 1.1alpha19

### The package

Get the database:

```

```

```
git clone https://gitlab.ebi.ac.uk/saary/eukcc
cd eukcc
python3 setup.py sdist bdist_wheel
pip3 install dist/eukcc-0.0.1.tar.gz
```

## Run it

```
EukCC.py --help
```

### Example

Run EukCC on a bunch of bins:

```
EukCC.py -c config/dir outdir bins/*.fa
```



## FAQ

**GeneMark-ES throws a lot of errors**
Yes, id does. If it does not start at all try the 32 bit version 
(with the 32 bit version licence). This worked for me.




