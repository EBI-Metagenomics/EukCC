=====================
Bin merging
=====================

If EukCC is run in ``folder`` mode, it can try to merge
two or more bins to create a refined/merged version
of increased completeness.

For this you can and should pass paired read information to EukCC, so
only bins linked by at least 100 reads (``--min_links``) are considered
for merging. This greatly improves speed and accuracy.

Preparing your linked reads
=============================
If you have paired-end read data you should create a sorted alignment. 
If you have multiple read files, you can create multiple BAM files.

For this you will need the contigs that were used to create these bins.
Alternatively merge all bins into a pseudo-assembly file.

.. code-block:: shell

    cat binfolder/*.fa > pseudo_contigs.fasta
    bwa index pseudo_contigs.fasta
    bwa mem -t 8 pseudo_contigs.fasta reads_1.fastq.gz reads_2.fastq.gz  | 
        samtools view -q 20 -Sb - | 
        samtools sort -@ 8 -O bam - -o alignment.bam
    samtools index alignment.bam

You can then create a bin_linking table by using the EukCC provided script:

.. code-block:: shell

    binlinks.py  --ANI 99 --within 1500 \
        --out linktable.csv binfolder alignment.bam


If you have multiple bam files, pass all of them to the script (e.g. `*.bam`).

You will obtain a three column file (bin_1,bin_2,links).

Merging bins
=====================================
You can then launch EukCC on the same binfolder like so:


.. code-block:: shell

    eukcc folder \
        --out outfolder \
        --threads 8  \
        --links linktable.csv \
        binfolder

EukCC will first run on all bins individually. Bins with at least 50%
completeness are then combined with smaller bins (less than 50%) that are
linked to them by at least ``--min_links`` paired end reads (default: 100).
Without a link table, every small bin is tried.

A merge is kept if completeness goes up by at least ``--improve_percent``
(default: 10) and by more than ``--improve_ratio`` (default: 5) times the
increase in contamination.

Merged bins can be found in the output folder.

To only estimate the quality of each bin, without merging, use
``--no_merge``:

.. code-block:: shell

    eukcc folder --out outfolder --threads 8 --no_merge binfolder

.. warning::
    Meging more than two bins. So setting ``--n_combine`` to anything above 1 is experimental and not yet recommended. We had very good results with merging two bins.

Example Dataset
=======================================
I created example data to test this based on the Lichen Study ERP123954. Bins were created using CONCOCT but any binner 
with no prokaryotic bias works.

.. code-block:: shell

    wget ...
    gunzip eukcc_example_folder_GT57.zip

    # Use at least a couple threads to speed it up
    eukcc folder --threads 6 \
        --out output \
        --links eukcc_example_folder_GT57/links.csv \
        --n_combine 1 \
        eukcc_example_folder_GT57/bins


