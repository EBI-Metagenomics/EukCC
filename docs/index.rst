.. EukCC documentation master file, created by
   sphinx-quickstart on Thu Oct 10 10:56:39 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to EukCC's documentation!
=================================

EukCC is a command line too written in python3 to estimate completeness and contamination
of novel eukaryotic MAGs.

.. If you use it please cite:

The project is hosted on github: https://github.com/openpaul/EukCC


Installation
---------------

Install dependencies
~~~~~~~~~~~~~~~~~~~~~~
EukCC depends on hmmer, Genemark-ES and pplacer. So you will need to install 
them before launching EukCC. It is best to create a conda enviorement
and then install Genemark-ES manually.


**With conda**

.. code-block:: shell

   conda create -n eukcc -c bioconda -c biocore hmmer pplacer python=3.7
   # we recommend installing ete3 from conda:
   conda activate eukcc
   conda install -c etetoolkit ete3 

To install GeneMark-ES you will need to install perl and certain perl packages:

.. code-block:: shell

   conda activate eukcc
   conda install -c anaconda perl \
                 -c bioconda perl-app-cpanminus 
   cpanm YAML
   cpanm Logger::Simple
   cpanm Parallel::ForkManager.pm

Once you installed these dependencies you need to download GeneMark-ES
from http://exon.gatech.edu/GeneMark/license_download.cgi and deposit the
license key in your home directory as `.gm_key`.

Make sure `gmes_petap.pl` is availiable from your $PATH:

.. code-block:: shell

   which gmes_petap.pl
   #~/software/genemark_es/gmes_petap.pl



Install EukCC
~~~~~~~~~~~~~~~~

**Install from github (access required)**

Fetch from github and install in the conda environment:

.. code-block:: shell
    
    git clone https://github.com/openpaul/eukcc
    cd eukcc
    conda activate eukcc
    pip3 install . 

EukCC requires python 3.7

**Get database from the EBI cluster**

This will only work for EBI users and this section will be removed in 
the future.

.. code-block:: shell
    
    cp /hps/nobackup2/production/metagenomics/saary/databases/eukcc/eukcc_db_20191023.tar.gz .
    tar -xzvf eukcc_db_20191023.tar.gz
    mv eukcc_db_20191023 eukccdb

Now that the database is fetched you can invoke eukcc with:

.. code-block:: shell
    
    eukcc --db eukccb MAG.fa


Quickstart and testing
~~~~~~~~~~~~~~~~~~~~~~

Assuming you installed eukcc and fetched the database you can start eukcc
and test for function with this command:

.. code-block:: shell

    wget -O testgenome.fa ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/251/995/GCF_002251995.1_ASM225199v2/GCF_002251995.1_ASM225199v2_genomic.fna.gz
    eukcc --db eukccdb  \
        --ncores 4 \
        --ncorespplacer 1 \
        --outdir eukcc_testgenome \
        testgenome.fa


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   index.rst
   tutorial.rst
   FAQ.rst



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
