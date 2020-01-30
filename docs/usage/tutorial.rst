=====================================================
Tutorial: Assemble, bin and asses eukaryotic bins
=====================================================

We assume for this tutorial, that you have sequenced 
and assembled your metagenomic reads without prior filtering
using your assembler of choice (metaSPAdes, flye, megahit, etc.).


Binning 
---------------------

To bin eukaryotic metagenomic sequences you need to use a binner with
no prior towards prokariotic sequences. Thus you will not be able
to use neither maxbin2 not metabat. We suggest the use of CONCOCT, but 
any unbiased binner will work.

Follow the github page of CONCOCT if you want:

https://github.com/BinPro/CONCOCT

We assume for this tutorial that your bins are in the folder 
fasta_bins
and end with .fa.
   
Optional: Filtering for euks
----------------------------------
After binning you may want to reduce the number of bins by looking for bins 
with large amounts of eukaryotic DNA. This step will bias your analysis
towards what will be detected by EukRep as such, but can be helpful
if you deal with a large number of bins.

You could for example filter for bins with a minimum number of 
eukaryotic bases and/or a ratio of eukaryotic to unclassified or bacterial
DNA. We noticed a ratio of at least 20 % eukaryotic DNA will already make sure
to exclude any bacterial bins.

This can be done using our helper script "filter_euk_bins.py"

Make sure you have EukRep (https://github.com/patrickwest/EukRep) installed.

.. code-block:: shell
   
    filter_euk_bins.py concoct_bins/*.fa 
    # optionally remove tmp folder again
    # rm -r tmp

This will create a file `assignment.csv` containing comma seperated 
columns with stats about the bin including one column 'eukaryotic', which
can help you decide.

The assignent file looks like:

.. csv-table:: assignment.csv
    :file: ../_static/assignment.csv
    :header-rows: 1


Try adjusting the ratios using the flags, if you are unhappy with the defaults.

.. code-block:: shell

   filter_euk_bins.py concoct_bins/*.fa --minbpeuks 500000 --eukratio 0.2


All options are described in the help:

.. code-block:: shell

   filter_euk_bins.py -h



Estimate bin completeness and contamination using EukCC
-------------------------------------------------------
Now we can launch EukCC for each eukaryotic bin. We suggest submitting 
EukCC to a compute cluster, in this example we use LSF bsub:

We assume that you located the eukccdb in your home folder, you might need
to adjust that path.

.. code-block:: shell

    while IFS=, read -r binpath binname passed bp_eukaryote bp_prokaryote bp_unassigned bp_sum
    do
        if [[ $passed == "True" ]]; then
            NAME=$binname
            bsub -M 40000 -J eukcc_$bin -n 16 \
                "eukcc --db $HOME/eukccdb --ncores 16 \
                 --plot \
                 --ncorespplacer 1 --outdir eukcc/$NAME $binpath"
        fi
    done < assignment.csv

In this example we allocate 40 Gb in memory for EukCC and use 16 cores for 
gene prediction and annotation. To reduce memory consumption we only use a 
single thread for pplacer.


EukCC ouput
#################

EukCC will create a structure like this:

.. code-block:: shell

    $ tree -L 2 eukcc/1.fa
        eukcc/1.fa
        ├── eukcc.tsv
        └── workfiles
            ├── gmes
            ├── hmmer
            └── pplacer

The main output is the file `eukcc.tsv`. It will contain predictions with up
to three sets chosen to best encompass the phylogenetic location of the bin.

If EukCC fails to process the bin, it is likely due to GeneMark-ES not being able
to predict any proteins. In this case you might consider predicting proteins ahead
of time using AUGUSTUS or another solution.
But often this means, the bin is of low quality and can be ignored when searching for high
quality MAGs.


.. csv-table:: eukcc.csv
    :file: ../_static/eukcc.csv
    :header-rows: 1

In this table up to three quality estimates are given. These are the three 
best sets found to estimate the quality of this MAG. The best set is the 
first row. It can be interesting to look at more than one set as sometimes 
a lower set gives a more robust estimate.

Most columns of the table will be self explanatory. We want to highlight a few 
critical ones:

- **n**: This is the number of protein-profiles used to estimate completeness
  and contamination

- **node**: This is the location in the reference tree and can be used to 
  see if different mags are located in the same area of the tree.

- **ngenomes**: This is the number of reference genomes used to construct the
  set used to estimate bin quality. A low number will suggest a less stable
  estimate.

- **nPlacements and cover**: `nPlacements` is the number of proteins placed 
  in the reference tree and `cover` is the number of these placed below the
  `node`-set used to estimate completeness.




