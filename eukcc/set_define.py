import os
import logging
from eukcc.eukcc import eukcc, eukcc_state, check_set_size


def define_set_w_genomes(args):
    state = eukcc_state(workdir=None, options=vars(args))
    # if db is not set, we check for env variable
    if state["db"] is None:
        if os.environ.get("EUKCC2_DB") is not None:
            state["db"] = os.environ.get("EUKCC2_DB")
            logging.debug(
                f"Defined db via env variable EUKCC2_DB as '{state['db']}'"
            )
        else:
            logging.error("No database was provided via --db or EUKCC2_DB env variable")
            exit(1)

    check_set_size(state)

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
        logging.warning("No marker gene set could be defined")
        E.terminate()
    E.state.save_state(filename=state["name"])
