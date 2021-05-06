from eukcc.eukcc import eukcc, eukcc_state
from eukcc.file import file
from eukcc.fasta import merge_fasta
import os


class bin:
    def __init__(self, state, workdir, fasta, protein=False):
        self.state = eukcc_state(workdir, state)
        self.state["fasta"] = fasta
        if protein:
            self.state["seqtype"] = "AA"
            self.state["predicted_proteins"] = "metaeuk"
        self.state["workdir"] = os.path.abspath(workdir)
        self.name = os.path.basename(fasta)
        state["name"] = self.name
        # adding caching so reruns are extra fast
        svp = os.path.join(self.state["workdir"], "save.json.gz")
        if os.path.exists(svp):
            self.state.load_state(svp)
        if self.state["estimated"] is not True:
            self.run_eukcc()
            self.state["estimated"] = True
            self.state.save_state(svp)

    def run_eukcc(self):
        """
        Run EukCC2 on this bin
        """
        # launch new EukCC instance
        E = eukcc(self.state)
        if E.state["seqtype"] == "DNA":
            # predict proteins
            E.predict_protein()
        else:
            # copy fna to faa state entry
            E.state["faa"] = E.state["fasta"]

        if E.placement() is None:
            return None
        # decide which db to use
        clade = E.determine_subdb()
        if clade != "base":
            E.state["clade"] = clade
            E.state["dbinfo"] = E.load_db(self.state["db"], clade=clade)
        # now that we loaded the new DB we can continue using the rest of the algorythm

        # pick marker set and placement in one step
        if E.pick_marker_set() is None:
            return None

        E.state["scmg_data"] = E.hmmsearch_scmg(
            E.state["workdir"],
            E.state["faa"],
            E.state["marker_set"]["profiles"],
            keep_hmm=True,
        )
        E.compute_quality(E.state["scmg_data"], E.state["marker_set"]["profiles"])

    def __str__(self):
        return "Bin: {}".format(self.name)


def merge_bins(parent, children, workdir):
    """
    Function to merge a sinlge parent with n children of class bin
    Returns an object of eukcc_state
    """
    # create an all_bins container
    all_bins = [parent]
    for c in children:
        all_bins.append(c)

    # define a name for this merged bin, this is used in
    # file paths and later for the output
    names = [x.name for x in all_bins]
    name = "merged" + "_".join(names)
    wd = os.path.join(workdir, name)

    # We use the childerns faa file to search for missing parent markers
    faa_path = os.path.abspath(os.path.join(wd, "{}.faa".format(name)))
    if not os.path.exists(faa_path):
        faas = [x.state["faa"] for x in children]
        file.isdir(wd)
        merge_fasta(faas, faa_path, seperator=None)

    # initial a new eukcc state usign the parent as a starting point
    state = eukcc_state(workdir, parent.state)
    state["faa"] = faa_path
    state["name"] = name
    state["workdir"] = os.path.abspath(wd)
    # Caching is done for relaunchs
    svp = os.path.join(state["workdir"], "save.json.gz")
    if os.path.exists(svp):
        state.load_state(svp)
        return state
    else:
        # Initialize a eukcc class
        E = eukcc(state)
        # hmmersearc is done using parents pressed hmms
        E.state["scmg_data"] = E.hmmsearch_scmg(
            E.state["workdir"],
            E.state["faa"],
            E.state["marker_set"]["profiles"],
            use_hmm=parent.state["estimate_hmm_path"],
        )
        # merge parent results with kids
        E.state["scmg_data"].extend(parent.state["scmg_data"])
        E.compute_quality(E.state["scmg_data"], E.state["marker_set"]["profiles"])
        # saving this computation for a quick reload
        E.state.save_state(svp)
        return E.state
