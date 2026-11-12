EukCC options explained
============================

EukCC on a genome (default)
--------------------------------
By default ``eukcc single`` assumes a nucleotide fasta file (genome or MAG)
and predicts its proteins with MetaEuk, using the MetaEuk database that
comes with the EukCC database.

.. code-block:: shell

    eukcc single --db eukccdb --out outfolder genome.fna


EukCC on proteins
---------------------------
If the proteins of a genome or MAG were already predicted with another
pipeline, pass the protein fasta file instead. EukCC detects the sequence
type on its own; use ``--AA`` (or ``--DNA``) to set it explicitly.

.. code-block:: shell

    eukcc single --db eukccdb --out outfolder --AA proteins.faa

Proteins passed this way are used as they are. Before submitting them,
remove any stop codons (``*``) and avoid spaces and the characters
``|``, ``[`` and ``]`` in the sequence names, as they cause issues with
the placement in the reference tree. EukCC only cleans the names of the
proteins it predicts with MetaEuk.


Choosing the database clade
---------------------------
Besides the base database, the EukCC database has clade specific
databases (``fungi``, ``protozoa`` and ``plant``). By default EukCC places
the genome in the base tree and picks the clade from the placement: each
placement votes for the clade database of its closest references, and the
genome goes to the clade with at least 60% of the votes, otherwise it stays in
the base database. If the clade of the genome is known, it can be set with
``--clade``:

.. code-block:: shell

    eukcc single --db eukccdb --out outfolder --clade fungi genome.fna

With ``--simple`` EukCC only uses the base database, which is faster but
not specific.


Marker gene sets
---------------------------
The marker gene set is chosen from the placement of the genome. These
options change which markers are accepted:

- ``--marker_prevalence``: percentage of the reference species of a set in
  which a marker must be found as single copy (default: 95).
- ``--set_size``: minimal number of markers in a set, smaller sets are not
  used (default: 20).
- ``--max_set_size``: maximal number of markers in a set (default: 500).
  Larger sets are reduced to a fixed random sample of this size, so runs are
  reproducible; 0 uses all markers.
- ``--set_number_species``: minimal number of reference species a set is
  defined on (default: 3). Reduce it if no set can be found.

To use exactly the same marker gene set for another genome, pass the
``savestate.json.gz`` written by a previous ``eukcc single`` run:

.. code-block:: shell

    eukcc single --db eukccdb --out outfolder \
        --use_placement previous_run/savestate.json.gz \
        genome.fna

A marker gene set can also be defined from a group of genomes with
``eukcc define_set``. It writes ``<name>.json.gz`` to the output folder,
which can then be passed to ``--use_placement``:

.. code-block:: shell

    eukcc define_set --db eukccdb --out sets --name myset genome_1.fna genome_2.fna


Other options
---------------------------
- ``--keep``: keep the working directory after the run.
- ``--extra``: also write ``scmg_marker_table.csv.gz``, with the hits for
  every marker of the set, searched without the gathering cutoffs.
- ``eukcc ncbi_update``: update the NCBI taxonomy database used by ete3.
- ``eukcc download_db``: download and extract the EukCC database into
  ``--out``, from the EBI FTP or from ``--url``.


EukCC -h
------------------------

.. literalinclude:: ../_static/eukcc-help.txt
    :language: text
