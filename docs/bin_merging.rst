=====================
Bin merging
=====================

If EukCC is run in ``folder`` mode, it can try to merge
two more more bins to create a refined/merged version
of increased completeness.


For this you can and should pass paired read information to EukCC. So
only bins linked by at least 100 (default) reads are considered
for merging. This greatly improves speed and accuracy.

Preparing your linked reads
=============================
If you have paired-end read data you should create a sorted alignment.

For this you will need the contigs that were used to create this bins.
Alternatively merge all bins into a pseudo-assembly file.

.. code-block:: shell

    cat binfolder/*.fa > pseudo_contigs.fasta
    bwa-mem -t 8 pseudo_contigs.fasta reads_1.fastq.gz reads_2.fastq.gz  | 
        samtools view -q 20 -Sb - | 
        samtools sort -@ 8 -O bam - -o alignment.bam
    samtools index alignment.bam

You can then create a bin_linking table by using the EukCC provided script:

.. code-block:: shell

    binlinks.py  --ANI 99 --within 1500 \
        --out linktable.csv binfolder alignment.bam


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

EukCC will fist run on all bins individually. It will then identify medium quality bins that are at least 50% compelete but not yet more than
100-``improve_percent``. 
It will then identify bins that are linked by at least 100 paired end reads to these medium quality bins. If after 
merging the quality score goes up this bin will be merged. 

Merged bins can be found in the output folder.

.. warning::
    Meging more than two bins. So setting ``--n_combine`` to anything above 1 is experimental and not yet recommended. We had very good results with merging two bins.
