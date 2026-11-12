==============================
Quickstart
==============================

Setup
=================================
EukCC is available to install from many sources. We recommend using the
BioContainers image (https://quay.io/repository/biocontainers/eukcc). These
images have no ``latest`` tag, so pick a version tag from the repository:

.. code-block:: shell

    docker pull quay.io/biocontainers/eukcc:<tag>
    # or
    singularity pull docker://quay.io/biocontainers/eukcc:<tag>

Other options are: 

**Bioconda** (https://anaconda.org/bioconda/eukcc)

.. code-block:: shell

    conda install -c conda-forge -c bioconda "eukcc>=2"

**pip**

.. code-block:: shell

    pip install eukcc

With pip, the external tools (MetaEuk, HMMER, pplacer and epa-ng) must be
installed separately.

You can also fetch the source code from GitHub: https://github.com/EBI-Metagenomics/EukCC


Database setup
------------------------------------------------
You will need to fetch the database for EukCC once. 

.. note::

    The database from version 1 is not compatible. So after updating EukCC make sure to update the database.

.. warning::

    Databases from ``eukcc2_db_ver_2026-11`` on are recommended with EukCC 2.3 or later.

    - **Opisthokont protists:** these databases keep the references of opisthokont protists
      (Choanoflagellata, Ichthyosporea, Filasterea, Rotosphaerida and *Corallochytrium*) in
      the ``protozoa`` database, and not in ``fungi``. EukCC 2.3 picks the clade database by a
      vote of the closest references of each placement, so these genomes go to ``protozoa``.
      Older versions send every genome with Opisthokonta in its lineage to ``fungi``, where
      these protists have no close references. With an older version, ``--clade protozoa``
      sends them to the right database.
    - **Mamiellophyceae:** EukCC 2.3 also roots the tree on a reference when it chooses the
      marker gene set. With older versions and the ``plant`` database of
      ``eukcc2_db_ver_2026-11``, some Mamiellophyceae genomes (*Bathycoccus*, *Ostreococcus*
      and *Micromonas*) are placed in streptophytes.

Databases available at http://ftp.ebi.ac.uk/pub/databases/metagenomics/eukcc/:

.. list-table::
    :header-rows: 1

    * - File
      - Last update
      - Size
      - EukCC version
    * - ``eukcc2_db_ver_2026-11.tar.gz``
      - 2026-11-01
      -
      - >= 2.3 (recommended)
    * - ``eukcc2_db_ver_1.2.tar.gz``
      - 2022-04-01
      - 6.0G
      - >= 2
    * - ``eukcc2_db_ver_1.1.tar.gz``
      - 2021-10-15
      - 5.9G
      - >= 2
    * - ``eukcc2_db_ver_1.tar.gz``
      - 2021-04-07
      - 5.9G
      - >= 2
    * - ``eukcc_db_v1.1.tar.gz``
      - 2020-01-27
      - 659M
      - 1

Fetching the database is as simple as:

.. code-block:: shell

    eukcc download_db --out eukccdb

``eukcc download_db`` downloads the database above and extracts it in the
``--out`` folder. Use ``--url`` for another version from the FTP. To do it by
hand instead:

.. code-block:: shell

    mkdir eukccdb
    cd eukccdb
    wget http://ftp.ebi.ac.uk/pub/databases/metagenomics/eukcc/eukcc2_db_ver_2026-11.tar.gz
    tar -xzvf eukcc2_db_ver_2026-11.tar.gz

.. note::

    ``eukcc2_db_tiny_fungi.tar.gz`` is a reduced database, built from 10 fungi genomes and
    300 PANTHER profiles. Its small size makes it suitable for tests, especially in CI/CD,
    but not for quality estimates. Run it on a fungi genome with ``--clade fungi``:

    .. code-block:: shell

        eukcc download_db --out eukccdb_tiny \
            --url https://ftp.ebi.ac.uk/pub/databases/metagenomics/eukcc/eukcc2_db_tiny_fungi.tar.gz
        eukcc single --db eukccdb_tiny/eukcc2_db_tiny_fungi --out outfolder --clade fungi genome.fna

If you want to forget about the database location you can add an ENV variable to your
`.bash.rc`:

.. code-block:: shell

    export EUKCC2_DB=/path/to/.../eukcc2_db_ver_2026-11


Running EukCC
===========================================

EukCC comes with two modes. You can run EukCC on a single bin or on a folder of bins.

Running EukCC on a single bin gives you the most options to tweak the 
parameters as you see fit. For most metagenomic workflows running EukCC on a folder
of bins might be the most simple thing to do.


EukCC on a single MAG
----------------------------------
We assume that you did set your ``$EUKCC2_DB`` to the correct location. 
If not please use the ``--db`` flag to pass the database to EukCC.

.. code-block:: shell

    eukcc single --out outfolder --threads 8 bin.fa 


EukCC will then run on 8 threads. You can pass nucleotide fastas
or proteomes to EukCC, plain or gzip compressed (``.gz``). It will
automatically try to detect if it has to predict proteins or not. 

By default it will never use more than a single thread for placing 
the genomes in the reference tree, to save memory. 


EukCC on a folder of bins
-----------------------------------------------

.. code-block:: shell

    eukcc folder --out outfolder --threads 8 bins

EukCC will use every file with the extension ``.fa``, ``.fna`` or ``.fasta``
in the folder as a bin, also gzip compressed (``.fa.gz``, ...). If any of these files is
not a valid FASTA file, empty files included, EukCC lists all of them and stops.

In folder mode EukCC will also try to refine bins automatically.
To learn more about that please see :doc:`bin_merging`.

For more options, such as protein input or choosing the database clade,
see :doc:`usage/options`.

