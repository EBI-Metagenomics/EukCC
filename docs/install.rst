
Installation
---------------

Install dependencies
~~~~~~~~~~~~~~~~~~~~~~
EukCC depends on hmmer, Genemark-ES and pplacer. So you will need to install 
them before launching EukCC. It is best to create a conda enviorement
and then install Genemark-ES manually.


**With conda**

.. code-block:: shell

   conda create -n eukcc -c bioconda -c biocore hmmer>3.2 pplacer python=3.7
   # we recommend installing ete3 from conda:
   conda activate eukcc
   conda install -c etetoolkit -c anaconda ete3 pyqt>5

To install GeneMark-ES you will need to install perl and certain perl packages:

.. code-block:: shell

   conda activate eukcc
   conda install -c anaconda perl \
                 -c bioconda perl-app-cpanminus 
   cpanm Test::Pod
   cpanm Logger::Simple
   cpanm Parallel::ForkManager.pm
   cpanm YAML
   cpanm Hash::Merge module

Once you installed these dependencies you need to download GeneMark-ES
from http://exon.gatech.edu/GeneMark/license_download.cgi and deposit the
license key in your home directory as `.gm_key`.

Make sure `gmes_petap.pl`  and `get_sequence_from_GTF.pl` are availiable from your $PATH:

.. code-block:: shell

   which gmes_petap.pl
   #~/software/genemark_es/gmes_petap.pl



Install EukCC
~~~~~~~~~~~~~~~~

**Install from github (access required)**

Fetch from github and install in the conda environment:

.. code-block:: shell
    
    git clone https://github.com/Finn-Lab/EukCC/ eukcc
    cd eukcc
    conda activate eukcc
    pip install .  # might be pip3 if you have two versions of python in your path

EukCC requires python 3.7

**Get database from the EBI cluster**

This will only work for EBI users and this section will be removed in 
the future.

.. code-block:: shell
    
    wget http://ftp.ebi.ac.uk/pub/databases/metagenomics/eukcc/eukcc_db_v1.1.tar.gz
    tar -xzvf eukcc_db_v1.1.tar.gz
    mv eukcc_db_20191023_1  eukccdb

Now that the database is fetched you can invoke eukcc with:

.. code-block:: shell
    
    eukcc --db eukccb MAG.fa


Quickstart and testing
~~~~~~~~~~~~~~~~~~~~~~

Assuming you installed eukcc and fetched the database you can start eukcc
and test for function with this command:

.. code-block:: shell

    wget -O testgenome.fa.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/251/995/GCF_002251995.1_ASM225199v2/GCF_002251995.1_ASM225199v2_genomic.fna.gz
    gunzip testgenome.fa.gz
    eukcc --db eukccdb  \
        --ncores 4 \
        --ncorespplacer 1 \
        --outdir eukcc_testgenome \
        testgenome.fa


