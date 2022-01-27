import os
import logging
from collections import Counter
from eukcc import eukcc
from eukcc.fasta import determine_type


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
    # validate fasta input
    if my_state["seqtype"] is None:
        my_state["seqtype"] = determine_type(my_state["fasta"])
        logging.info("Set sequence type to {}".format(my_state["seqtype"]))

    # initialise a new EukCC instance
    E = eukcc.eukcc(my_state)

    if E.state["seqtype"] == "DNA":
        # predict proteins
        E.predict_protein()
    else:
        # copy fna to faa state entry
        E.state["faa"] = E.state["fasta"]

    markers = E.search_markers(use_all=True)

    # reduce to singletons
    found = Counter([row["query"] for row in markers])
    single = set()
    for profile, number in found.items():
        if number == 1:
            single.add(profile)

    logging.info("Done processing genome")
    return single
