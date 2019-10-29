=====================================================
Tutorial: Assemble, bin and asses eukaryotic bins
=====================================================

We assume for this tutorial, that you have sequenced 
and assembled your metagenomic reads without prior filtering
using your assembler of choice (metaSPAdes, flye, megahit, etc.).

We assume the scaffolds to be in the folder.


Binning 
---------------------

To bin eukaryotic metagenomic sequences you need to use a binner with
no prior towards prokariotic sequences. Thus you will not be able
to use neither maxbin2 not metabat. We suggest the use of CONCOCT, but 
any unbiased binner will work.

.. code-block:: shell

   concoct .. bin 
   
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

.. code-block:: shell
   
   filter_euk_bins.py fasta_bins/*.fa

This will create a file `asignment.csv` containing comma seperated 
columns with stats about the bin including one column 'eukaryotic', which
can help you decide.

The assignent file looks like:

   path,binname,eukaryotic,eukbases,bacbases,unassigned,sum
   fasta_bins/0.fa,0.fa,True,2593713,0,0,2593713
   fasta_bins/1.fa,1.fa,False,36682,394058,3639,434379

Try adjusting the ratios using the flags, if you are unhappy with the defaults.

.. code-block:: shell

   filter_euk_bins.py fasta_bins/*.fa --eukratio 0.5 --bacratio 0.4


All options are described in the help:

.. code-block:: shell

   filter_euk_bins.py -h




Estimate bin completeness and contamination using EukCC
-------------------------------------------------------


.. code-block:: shell

   eukcc --db $HOME/eukccdb/ -n 8 --outdir eukcc/bin.0/ fasta_bins/bin.0.fa

Optimizing memory
##################
Running EukCC will require serveral Gb of memory (>50 Gb). This is due to
the usage of pplacer. 

In a first step, it will run GeneMark-ES, which depending on the 
size of the bin and the number of cores provided can take up to several
hours and uses very little memory. 

If you want to streamline your pipeline, you can use EukCC to only run 
GeneMark-ES first (on a low memory machine) and then relaunch with the same 
parameters on a higher memory machine. 
