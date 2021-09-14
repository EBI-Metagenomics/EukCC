import os
import logging
from collections import Counter
from eukcc import eukcc


def singletons_from_genome(fasta, state):
    """
    search a genome for all markers in a clade
    """
    logging.info("Running for {}".format(fasta))

    workdir = os.path.join(
        state["workdir"], "genome_annotation", os.path.basename(fasta)
    )
    my_state = eukcc.eukcc_state(workdir, state)
    my_state["fasta"] = fasta
    E = eukcc.eukcc(my_state)
    E.predict_protein()
    markers = E.search_markers(use_all=True)

    # reduce to singletons
    found = Counter([row["query"] for row in markers])
    single = set()
    for profile, number in found.items():
        if number == 1:
            single.add(profile)

    logging.info("Done processing genome")
    return single
