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


Quickstart
~~~~~~~~~~~~~~~~~~~~~~

The database is currently located under:

.. code-block:: shell

    # 14. october 2019
    /hps/nobackup2/production/metagenomics/saary/databases/eukcc/balanced5_einsi


And then launch EukCC with:

.. code-block:: shell

   eukcc --db eukccdb MAG.fa


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   index.rst
   tutorial.rst



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
