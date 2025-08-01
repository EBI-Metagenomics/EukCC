==============================
Quickstart
==============================

Setup
=================================
EukCC is avaliable to install frommany soruces. We reccomend using our 
Docker container:

.. code-block:: shell

    docker pull microbiomeinformatics/eukcc

Other options are: 

**Bioconda** (https://anaconda.org/bioconda/eukcc)

.. code-block:: shell

    conda install -c conda-forge -c bioconda "eukcc>=2"

    pip install eukcc

You can also fetch the source code from GitHub: https://github.com/EBI-Metagenomics/EukCC


Database setup
------------------------------------------------
You will need to fetch the database for EukCC once. 

.. note::

    The database from version 1 is not compatible. So after upadting to EukCC make sure to upadte the database.

Fetching the database is as simple as:

.. code-block:: shell

    mkdir eukccdb
    cd eukccdb
    wget http://ftp.ebi.ac.uk/pub/databases/metagenomics/eukcc/eukcc2_db_ver_1.2.tar.gz
    tar -xzvf eukcc2_db_ver_1.2.tar.gz

If you want to forget about the database location you can add an ENV variable to your
`.bash.rc`:

.. code-block:: shell

    export EUKCC2_DB=/path/to/.../eukcc2_db_ver_1.2


Running EukCC
===========================================

EukCC comes with two modes. You can run EukCC on a single bin or on a folder of bins.

Running EukCC on a single bin gives you the most options to tweak the 
parameters as you see fit. For most metagenomic workflows running EukCC on a folder
of bins might be the most simple thing to do.


EukCC on a single MAG
----------------------------------
We assume that you did set you ``$EUKCC2_DB`` to the correct location. 
If not please use the ``--db`` flag to pass the database to EukCC.

.. code-block:: shell

    eukcc single --out outfolder --threads 8 bin.fa 


EukCC will then run on 8 threads. You can pass nucleotide fastas
or proteomes to EukCC. It will automatically try to detect if it
has to predict proteins or not. 

By default it will never use more than a single threads for placing 
the genomes in the reference tree, to save memory. 


EukCC on a folder of bins
-----------------------------------------------

.. code-block:: shell

    eukcc folder --out outfolder --threads 8 bins

EukCC will assume that the folder contains files with the suffix ``.fa``. If 
that is not the case please adjust the parameter.

In folder mode EukCC will also try to refine bins automatically. 
To learn more about that please see 

