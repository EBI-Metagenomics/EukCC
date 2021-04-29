
Installation
---------------

EukCC is avaliable to install from bioconda (https://anaconda.org/bioconda/eukcc)
or pypi (https://pypi.org/project/eukcc/).

You can also fetch the source code from GitHub: https://github.com/Finn-Lab/EukCC/

Install via conda
~~~~~~~~~~~~~~~~~~~~~~
Installing EukCC via conda will install all dependencies 
**except for GeneMark-ES**.  If you want to use GeneMark-ES
you will need to download a 64 bit version with the license here:
http://exon.gatech.edu/GeneMark/license_download.cgi

.. code-block:: shell

    conda install -c bioconda  -c conda-forge eukcc perl-app-cpanminus
    # now install perl dependencies for GeneMark-ES using cpanminus
    env PERL5LIB="" PERL_LOCAL_LIB_ROOT="" PERL_MM_OPT="" PERL_MB_OPT="" $CONDA_PREFIX/bin/cpanm inc::Module::Install::DSL Hash::Merge MCE::Mutex FindBin Test::Pod Logger::Simple  Parallel::ForkManager.pm YAML Math::Utils

Once you installed these dependencies you need to download GeneMark-ES
from http://exon.gatech.edu/GeneMark/license_download.cgi and deposit the
license key in your home directory as `.gm_key`.

.. code-block:: shell

    tar -xzvf gmes_linux_64.tar.gz
    cp gmes_linux_64/* ~/local/bin/
    
    # extract the key file
    zcat gm_key.gz > ~/.gm_key

Make sure `gmes_petap.pl`  and `get_sequence_from_GTF.pl` are availiable from your $PATH:

.. code-block:: shell

   which gmes_petap.pl
   #~/software/genemark_es/gmes_petap.pl


**Get database from the EBI cluster**

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

