import os
import logging
from eukcc.eukcc import eukcc, eukcc_state


def define_set_w_genomes(args):
    state = eukcc_state(workdir=None, options=vars(args))
    # if db is not set, we check for env variable
    if state["db"] is None:
        if os.environ.get("EUKCC2_DB") is not None:
            state["db"] = os.environ.get("EUKCC2_DB")
            logging.debug(
                "Defined db via env variable EUKCC2_DB as '{}'".format(state["db"])
            )
        else:
            logging.error("No database was provided via --db or EUKCC2_DB env variable")
            exit(202)

    if state["max_set_size"] < state["set_size"]:
        logging.warning("Max set size cant be smaller than min set size")
        exit(222)

    if state["seqtype"] is None:
        logging.debug("Setting seqencetype to DNA, use --AA to pass proteomes")
        state["seqtype"] = "DNA"
    state["comparison_genomes"] = state["genomes"]
    state["ignore_tree"] = True
    state["fasta"] = state["genomes"][0]
    # use base as workdir folder
    state["workdir"] = os.path.join(state["out"], "workdir")
    # launch new EukCC instance
    E = eukcc(state)
    if E.pick_marker_set() is None:
        E.terminate(1)
    E.state.save_state(filename=state["name"])
