EukCC options explained
============================

EukCC in genome (default) mode
--------------------------------
When launching EukCC without special 
parameters, it will assume that a 
genomic fasta file was passed as input.
Thus GeneMark-ES will be used 
to predict proteins. 

.. code-block:: shell

    eukcc --db eukccdb -o . genome.fna


EukCC in protein mode
---------------------------
If proteins for a genome or MAG were already predicted
using a nother pipeline, EukCC can be used to 
estimate the completeness and contamination.

For this EukCC requires at least the proteins as Fasta
file. Internally EukCC ignores repeated proteins that
occure very close to each other on a genomic level. This
is due to common gene prediction errors and subsequent 
duplicated hits with profile hmms. Thus it is possible
to pass the genomic coordinates for each protein as a bed
file to EukCC. Such a file can easily be prepared from 
gtf or gff files. If no bed file is provided, this step 
will be skipped.

The bed file needs to be a simple five column bed file containing the contig, start and end of the gene, 
the strand and finally the name of the gene. Exons and introns are ignored as only start and end of the coding regions matter. 

Befor submissint predicted proteins make sure to remove any stop codons (*), as they cause issues with pplacer (https://github.com/Finn-Lab/EukCC/issues/14).


.. code-block:: shell

    eukcc --db eukccdb -o . \
          --protein genome.faa \
          --bed coordinates.bed
        


EukCC using pygmes
---------------------------
GeneMark-ES uses a self training step to generate a model
for the provided genome. Somtimes highly fragmented or
incomplete genomes can fail to provided enough training data
for GeneMark-ES model creation to succeed. 

In these cases EukCC can rely on pygmes to select a suitable 
model from previous analysed genomes. This allows to estimate 
genome completeness also for very incomplete genomes.

For pygmes to select a suitable model the user needs to provide
a diamond data base with taxonomic information.
Such a database can be downloaded here:

.. code-block:: shell

    wget -O uniref50_pygmes.dmnd http://ftp.ebi.ac.uk/pub/databases/metagenomics/eukcc/uniref50_20200213_tax.dmnd

EukCC can then be launched with:


.. code-block:: shell

    eukcc --db eukccdb -o . \
        --pygmes \
        --diamond uniref50_pygmes.dmnd \
        genome.fna

    
EukCC -h
------------------------


.. include:: ../_static/eukcc-help.txt

