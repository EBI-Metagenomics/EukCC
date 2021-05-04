import logging
import os
import shutil
import json
from collections import defaultdict, Counter
import hashlib
import jsonpickle
import csv
import gzip
from eukcc.exec import (
    hmmsearch,
    hmmalign,
    epa_split,
    epa_ng,
    guppy,
    hmmpress,
    hmmfetch,
    metaeuk,
)
from eukcc.base import (
    which,
    read_hmmer,
    horizontal_concatenate,
    load_tax_info,
    load_SCMGs,
    compute_silent_contamination,
    hmmer_cutoffs,
)
import eukcc.fasta as Fasta
from eukcc.fasta import check_alignment
from eukcc.file import file
from eukcc.treehandler import tree_sets, hard_set_computation, tax_LCA
from glob import glob


default_vals = {
    "clade": "base",
    "seqtype": "DNA",
    "set_number_species": 3,
    "marker_prevalence": 98,
    "set_size": 20,
    "max_set_size": 500,
    "set_selection": "best_guess",
}


class eukcc_state:
    """
    This is a class to maintain the state of eukcc
    This allows us to write intermediate results to disk
    And thus resume again
    Usually you dont have to worry about this
    """

    def __init__(self, workdir="", options=None):
        """
        Only parameter that is required is the workdir,
        as that could be the folder with the intermediate results
        """
        self.opt = defaultdict(lambda: None)
        self.merge(default_vals)
        if options is not None:
            self.merge(options)
        self.opt["workdir"] = workdir

    def keys(self):
        return self.opt.keys()

    def __setitem__(self, key, value):
        self.opt[key] = value

    def __getitem__(self, key):
        return self.opt[key]

    def get(self, key):
        if key in self.opt.keys():
            return self.opt[key]
        else:
            return None

    def __str__(self):
        s = ""
        for key, value in self.opt.items():
            s += "{}: {}\n".format(key, value)
        return s

    def checkpoint(self, value):
        self.chkpnt = value
        self.opt["last_checkpoint"] = value

    def save_state(self, path=None, filename="savestate"):
        if path is None:
            jo = os.path.join(self.opt["out"], "{}.json.gz".format(filename))
        else:
            jo = path
        with gzip.open(jo, "wb") as fout:
            d = json.dumps(jsonpickle.encode(self.opt))
            fout.write(d.encode())
        return

    def load_state(self, path=None):
        if path is None:
            path = os.path.join(self.opt["out"], "savestate.json.gz")
        if not os.path.exists(path):
            logging.warning("Can not load save as file does not exist")
        else:
            with gzip.open(path, "r") as fin:
                data = json.loads(fin.read().decode("utf-8"))
                opt = jsonpickle.decode(data)
                self.merge(opt)
        return

    def merge(self, d):
        if type(d) is eukcc_state:
            opt = d.opt
        else:
            opt = d
        for key, value in opt.items():
            self.opt[key] = value


def make_absolute(path):
    return os.path.abspath(path)


class eukcc:
    """
    Evaluate completeness and contamination of eukaryotic
    genomes or MAGs
    """

    def __init__(self, state):
        # state is kept internally
        self.state = state
        self.state.checkpoint("launch")
        for key in ["db", "fasta", "out"]:
            if key in state.keys():
                self.state[key] = make_absolute(state[key])

        logging.debug("Using database located at: {}".format(self.state["db"]))
        if self.state["workdir"] is None or self.state["workdir"] == "":
            self.state["workdir"] = os.path.join(state["out"], "workdir")
        file.isdir(self.state["workdir"])

        # check if the database is accessible
        self.state["loaded_dbs"] = {}
        self.state["dbinfo"] = self.load_db(state["db"], "base")
        if state["clade"] != "base":
            self.state["dbinfo"] = self.load_db(state["db"], state["clade"])
        logging.debug(
            "Using database '{}' with PANTHER version {}".format(
                self.state["dbinfo"]["name"], state["dbinfo"]["panther_version"]
            )
        )

        # validate AA type
        self.state["seqtype"] = state["seqtype"].upper()
        if self.state["seqtype"] not in ["AA", "DNA"]:
            raise ValueError("seqtype should be either DNA or AA")

        if Fasta.validate_fasta(self.state["fasta"]) is False:
            raise ValueError("The provided fasta is corrupt or maybe gzipped?")

        if self.state["set_number_species"] < 3:
            logging.warning("You included less than 3 species in the set creation, this might lead to unstable results")
        if self.state["marker_prevalence"] < 95:
            logging.warning("You selected a marker prevalence below 95%, this can lead to wrong results")

        if self.state["set_size"] > state["max_set_size"]:
            self.state["max_set_size"] = state["set_minimal_size"]
            logging.warning(
                "Your minimal set size is bigger than the maximial set size. Thus we increased the maximal set size to fit the minimal set size. Both are now {}".format(
                    self.state["max_set_size"]
                )
            )

        logging.debug("Using {threads} threads".format(threads=self.state["threads"]))

        if self.state["rerun"]:
            logging.debug("deleting workdir to rerun all of the analysis")
            self.remove_workdir(self.state["workdir"])

        # check dependencies
        self.check_dependencies()

    def placement(self):
        # place in reference tree using epa, if we
        # do not get taxids to orient us
        if self.state["taxids"] is not None:
            logging.info("Skipping placement: Going to use taxids as seed for set construction")
            # infer placement from taxids
            self.state["placement"] = self.ncbi_place(self.state)
            self.state["placement"]["tog"] = self.state["dbinfo"]["files"]["backbone_tree"]
        else:
            logging.info("Searching for marker genes in {} database".format(self.state["clade"]))
            self.state["marker_genes"] = self.search_markers()
            if len(self.state["marker_genes"]) == 0:
                logging.info("No placement marker genes found.")
                return None
            logging.info(
                "Found {} marker genes, placing them in the tree using epa-ng".format(len(self.state["marker_genes"]))
            )
            pl = self.run_EPA()
            if pl is None:
                return None
            else:
                self.state["placement"] = pl
        return True

    def pick_marker_set(self):
        if self.state["use_placement"] is not None:
            logging.debug("Replying on provided placement file")
            old_state = eukcc_state(workdir=os.path.join(self.state["workdir"], "loading"))
            # load the old info
            old_state.load_state(self.state["use_placement"])
            self.state["marker_set"] = old_state["marker_set"]
            return self.state

        # figure placement
        if self.placement() is None:
            return None

        if self.state["ignore_tree"] and self.state["placement"]["genomes"] is not None:
            logging.info("Will ignore tree and compare to provided genomes directly")
            ms = hard_set_computation(
                self.state["dbinfo"]["files"]["scmgs"],
                genomes=self.state["placement"]["genomes"],
                prevalence=self.state["marker_prevalence"],
                atmost=self.state["max_set_size"],
                set_size=self.state["set_size"],
            )
            if ms is None:
                logging.error(
                    "No marker gene set could be found with these settings \U0001F641 \nChange your parameters and try again?"
                )
                return None
            else:
                self.state["marker_set"] = ms
        else:
            logging.info("Automatically locating best SCMG set")
            tree = tree_sets(
                self.state["placement"]["tog"],
                self.state["placement"],
                self.state["dbinfo"]["files"]["scmgs"],
                set_species=self.state["set_number_species"],
                set_size=self.state["set_size"],
                set_prevalence=self.state["marker_prevalence"],
                set_atmost=self.state["max_set_size"],
                dynamic_root=self.state["dynamic_root"],
                use_ncbi=self.state["use_ncbi"],
                taxinfo=self.state["dbinfo"]["files"]["taxinfo"],
                set_selection=self.state["set_selection"],
            )
            if tree.marker_set is None:
                logging.error(
                    "No marker gene set could be found with these settings \U0001F641 \nChange your parameters and try again?"
                )
                return None
            if len(tree.marker_set.covered) < (0.5 * len(tree.marker_set.all_places)):
                n = len(tree.marker_set.covered)
                m = len(tree.marker_set.all_places)
                logging.warning(
                    "\U00002757The choosen marker gene set is supported by only half ({}/{}) of the alignments. This generally is an unstable estimate.".format(
                        n, m
                    )
                )
            self.state["marker_set"] = tree.marker_set.to_dict()
            return True

    def predict_protein(self):
        # fasta_faa = self._pygmes(workdir, fasta, ncores=ncores)
        logging.debug("Fasta is of DNA type, proteins will be predicted")
        if self.state["use_gmes"] is True:
            logging.debug("Predicting proteins using GeneMark-ES")
            fasta_faa = self._pygmes(
                self.state["workdir"],
                self.tate["fasta"],
                diamond=self.state["dbinfo"]["files"]["diamond"],
                ncores=self.state["threads"],
            )
        else:
            logging.debug("Predicting proteins using MetaEuk")
            fasta_faa = self._metaeuk(
                self.state["workdir"],
                self.state["fasta"],
                self.state["loaded_dbs"]["base"]["files"]["metaeukdb"],  # only the base db has the metaeuk db
                ncores=self.state["threads"],
            )
        self.state["faa"] = fasta_faa

    def _metaeuk(self, workdir, fasta, metaeukdb, ncores=1):
        """
        Call metaeuk with the fasta while using a nice workdir
        """
        wd = os.path.join(workdir, "metaeuk")
        prefix = "{}_metaeuk".format(os.path.basename(fasta).rsplit(".", 1)[0]).replace(".", "_")
        logging.debug("Metaeuk prefix set as {}".format(prefix))
        faafile = "{}.fas".format(prefix)
        cleanfaafile = "{}_cleaned.faa".format(prefix)

        logging.debug("Launching metaeuk to produce {}".format(faafile))
        metaeuk(
            "metaeuk",
            wd,
            [fasta, metaeukdb],
            [faafile, cleanfaafile],
            prefix=prefix,
            cores=ncores,
        )

        self.state["predicted_proteins"] = "metaeuk"
        # remove all tmp files
        file.delete_but(wd, keep=[faafile, cleanfaafile])

        return os.path.join(wd, cleanfaafile)

    def _pygmes(self, workdir, fasta, diamond, ncores=1):
        """
        launch pygmes
        """
        from pygmes import pygmes

        wd = os.path.join(workdir, "pygmes")
        fasta_faa = os.path.join(wd, "predicted_proteins.faa")
        if file.isnewer(fasta, fasta_faa):
            pygmes(fasta, wd, db=diamond, clean=True, ncores=ncores)
        # delete all tmp files by default
        file.delete_but(wd, "predicted_proteins.faa")
        return fasta_faa

    def structure_result(self, **args):
        args["found_SCMGs"] = [x["query"] for x in args["found_SCMGs_raw"]]

        # make sure we only has lists, else json write will not work
        for key in args.keys():
            if type(args[key]) is set:
                args[key] = list(args[key])

        # store this object as json for now
        # so one could use it
        # in the future this logic might change
        # this is really open to debate and input
        # for now it works storing the json I assume
        # we store it as .gz because we can
        jo = os.path.join(args["outdir"], "results.json.gz")
        with gzip.open(jo, "wb") as fout:
            d = json.dumps(args)
            fout.write(d.encode())
        return args

    def determine_subdb(self):
        """
        assuming a placement was done, we try to figure out what tax clade we deal with
        """
        pl = self.state["placement"]
        tree = pl["tog"]
        places = [x["n"][0] for x in pl["placements"]]
        info = self.state["dbinfo"]["files"]["taxinfo"]
        lng = tax_LCA(tree, info, places)
        clade = "base"
        if len(lng) < 3:
            logging.warning("Could not determine a suffciently deep LCA, using base database")
        else:
            if len(set(["33154", "4751"]) & set(lng)) > 0:
                clade = "fungi"
            elif len(set(["33090"]) & set(lng)) > 0:
                clade = "plant"
            elif len(set(["protist_common"]) & set(lng)) > 0:
                clade = "protozoa"
        if clade != "base":
            logging.info("Genome belongs to clade: {} (Best TaxID: {})".format(clade, lng[-1]))
        else:
            logging.warning(
                "Could not narrow down the clade, using base database. Please provide clade if known via --clade {fungi,plant,protozoa}"
            )

        return clade

    def remove_workdir(self, workdir):
        """
        Remove folder recursivly. This is called if a rerun should be forced and if temoprary files are to be removed

        :param workdir: path to folder to remove
        """
        # no need to delete non existing folders
        if not os.path.exists(workdir):
            return
        logging.debug("Deleting workdir recursivly")
        if os.path.isdir(workdir):
            shutil.rmtree(workdir)
        elif os.path.exists(workdir):
            raise OSError("{} is not a folder but should be, please move or remove this file".format(workdir))

    def ncbi_place(self, state):
        """
        Given taxids we can find all reference species used to construct the
        backbone tree with overlapping taxonomy.
        Requires database to be loaded.

        :param taxids: List or set of taxids
        :return: dict with placements
        """
        # parse, so we allow comma and spaces
        taxa = []
        for t in self.state["taxids"]:
            if "," in t:
                taxa.extend(t.split(","))
            else:
                taxa.append(t)
        # convert taxa in a set of strings
        taxa = set([str(t) for t in taxa])
        info = load_tax_info(state["dbinfo"]["files"]["taxinfo"])
        # find all nodes that intersect with the taxids
        nodes = set()
        for node, lng in info.items():
            if len(taxa & set(lng)) > 0:
                nodes.add(node)
        # make sure these are also in our SCMG set
        scmgs = load_SCMGs(state["dbinfo"]["files"]["scmgs"])
        nodes = nodes & set(scmgs.keys())

        placements = [{"n": x} for x in nodes]
        logging.info("Located {} species corresponding to the provided taxids".format(len(nodes)))
        return {"placements": placements, "genomes": nodes}

    def compute_quality(self, data, profiles):
        """
        Compute completeness and contamination.

        :param data: read_hmmer data object with found hits in query column
        :param profiles: Expected SCMGs
        :return: dict of completenss and contamination
        """
        expected = set(profiles)
        found = Counter([x["query"] for x in data])
        found_marker = set(found.keys()) & expected
        if logging.DEBUG >= logging.root.level:
            missing_marker = set(expected - found.keys())
            logging.debug("Missing these expected markers: {}".format(missing_marker))
        multiton_marker = set([k for k, v in found.items() if v != 1]) & expected

        completeness = round(len(found_marker) / len(expected) * 100, 2)
        contamination = round(len(multiton_marker) / len(expected) * 100, 2)
        quality = {"completeness": completeness, "contamination": contamination}
        logging.info("Completeness: {}".format(completeness))
        logging.info("Contamination: {}".format(contamination))

        if self.state.get("predicted_proteins") == "metaeuk":
            _silent_contamination = compute_silent_contamination(self.state["fasta"], self.state["faa"], data)
            quality["silent_contamination"] = round(
                100 - 100 * _silent_contamination["bp_w_hits"] / _silent_contamination["total_bp"],
                2,
            )
            logging.info("Max silent contamination: {}".format(quality["silent_contamination"]))

        self.state["quality"] = quality
        self.state["quality"]["N50"] = Fasta.N50(self.state["fasta"])
        self.state["quality"]["N50_AA"] = Fasta.N50(self.state["faa"])

    def hmmsearch_scmg(self, workdir, fasta_faa, profiles, keep_hmm=False, use_hmm=None, cut_ga=True):
        """
        Fetch and compress and search markers specified by profiles list.
        Will call hmmfetch, hmmpress and hmmsearch internally

        :param workdir: folder to create tmp files under
        :param fasta_faa: predicted proteome
        :param profiles: list or set of profiles, must be in reference database
        """

        # creating a hash from the profiles
        # so set adjustments will lead to new results
        profiles = list(profiles)
        profiles.sort()
        m = hashlib.sha256()
        # concatenate all parameters to a single hash
        hash_str = "-".join(profiles)
        hash_str = hash_str + str(cut_ga)

        m.update(hash_str.encode())
        key = m.hexdigest()[0:12]

        logging.debug("Searching for selected markers")
        wd = os.path.join(workdir, "hmm", "singlecopy", self.state["clade"])
        file.isdir(wd)
        outfile = "found_markers_{}.tbl".format(key)
        cutoff_file = "found_markers_{}_GA.csv".format(key)

        # only run if outfile is older than fasta_faa
        if file.isnewer(fasta_faa, os.path.join(wd, outfile)) or file.isnewer(fasta_faa, os.path.join(wd, cutoff_file)):
            if use_hmm is None:
                # fetch and press hmms from database
                hmmfile = self.hmmfetch(wd, profiles)
                logging.debug("Pressing hmms")
                hmmpress("hmmpress", wd, hmmfile, "selected_hmms.hmm.h3m")
            else:
                hmmfile = use_hmm

            logging.info("Searching fasta for selected markers")
            hmmsearch(
                "hmmsearch",
                workdir=wd,
                infiles=[hmmfile, fasta_faa],
                outfiles=[outfile],
                cores=self.state["threads"],
                cut_ga=cut_ga,
            )

            # read GA values from fetched hmm file, needed for some stats in the end
            cutoff_file = hmmer_cutoffs(hmmfile, os.path.join(wd, cutoff_file))
            if keep_hmm:
                self.state["estimate_hmm_path"] = hmmfile
            else:
                # we can remove hmms and pressed hmms
                file.delete_but(
                    wd,
                    keep=[os.path.basename(x) for x in glob(os.path.join(wd, "found_markers_*"))],
                )

        return read_hmmer(os.path.join(wd, outfile), os.path.join(wd, cutoff_file))

    def hmmfetch(self, workdir, profiles):
        """
        Fetch profiles from database and copy them into a new
        reduced hmm file to be used for hmmscan

        :param workdir: Folder to create files in
        :param profiles: List or set of profile names to fetch
        """
        logging.debug("Fetching hmms")
        profile_file = os.path.join(workdir, "profiles.txt")
        with open(profile_file, "w") as fout:
            for p in profiles:
                fout.write("{}\n".format(p))
        outfile = "selected_hmms.hmm"
        hmmfetch(
            "hmmfetch",
            workdir,
            [self.state["dbinfo"]["files"]["hmm_db"], profile_file],
            outfile,
        )
        logging.debug("fetched hmms to file {}".format(outfile))
        return os.path.join(workdir, outfile)

    def terminate(self, code=0):
        """
        Terminate eukcc and write whatever we have to disk
        """
        exit(code)

    def guppy_tree(self, workdir, placement_file):
        """
        given a output file from epa-ng or pplacer will
        return a path to a file cotnaining a tree with the placed nodes in it

        :param wordir: temp folder
        :param placement_file: path to .jplace file
        :retrun: path to newick tree file
        """
        wd = os.path.join(workdir, "epa-ng", "guppy", self.state["clade"])
        guppy("guppy", wd, placement_file, "epa.tre")
        return os.path.join(wd, "epa.tre")

    def align_marker_genes(self, fasta, markers, workdir, output):
        """
        Align marker genes identified using hmmsearch to be passed on the EPA

        :params fasta: path to a protein fasta file
        :params markers: dictionary of marker genes obtained by read_hmmer
        :params workdir: Folder to write tmp files to
        :params output: Path to the final expected output file
        """
        wd = workdir
        file.isdir(wd)
        md = defaultdict(list)
        for row in markers:
            md[row["query"]].append(row["target"])
        # prots file
        prots = {}

        all_profiles = set(self.state["dbinfo"]["files"]["hmm_placement_place"].keys())
        for profile, ids in md.items():
            if profile not in self.state["dbinfo"]["files"]["hmm_placement_place"].keys():
                logging.error("Profile {} misses the hmm and alignment files".format(profile))
                exit(200)
            hmmfile = self.state["dbinfo"]["files"]["hmm_placement_place"][profile]["hmm"]
            alnfile = self.state["dbinfo"]["files"]["hmm_placement_place"][profile]["aln"]
            prot_file = os.path.join(wd, "{}.faa".format(profile))
            prot_file = Fasta.reduce_fasta(fasta, prot_file, ids)

            aligned = os.path.join(wd, "{}.fasta".format(profile))
            hmmalign(
                "hmmalign",
                wd,
                [alnfile, hmmfile, prot_file],
                [aligned],
            )
            # check that alignment has also non gap sites, else epa-ng will fail
            if check_alignment(aligned, prot_file):
                prots[profile] = aligned
                all_profiles.remove(profile)
            else:
                logging.debug("Skipped alignment, as it is only gaps")

        if len(prots.keys()) == 0:
            logging.debug("No alignments could be done")
            return None

        for profile in all_profiles:
            prots[profile] = self.state["dbinfo"]["files"]["hmm_placement_place"][profile]["aln"]

        logging.debug("Concatenating alignments")
        return horizontal_concatenate(output, [v for k, v in prots.items()], list(prots.keys()))

    def run_EPA(self):
        wd = os.path.join(self.state["workdir"], "epa-ng", self.state["clade"])
        file.isdir(wd)
        aligned_file = os.path.join(wd, "alignment.fasta")
        if file.isnewer(self.state["faa"], aligned_file):
            align_wd = os.path.join(wd, "align")
            aligned_file = self.align_marker_genes(
                self.state["faa"], self.state["marker_genes"], align_wd, aligned_file
            )
            if aligned_file is None:
                logging.debug("No marker genes aligned")
                return None
            # once the horizontal alignment is created we can delete
            # all small alignments
            file.remove(align_wd)

        # then split into query
        RS = epa_split(
            "epa-ng",
            wd,
            [self.state["dbinfo"]["files"]["backbone_alignment"], aligned_file],
            ["query.fasta", "reference.fasta"],
        )
        # Return False if placement failed
        if RS.success is False:
            return None

        # then epa place them
        R = epa_ng(
            "epa-ng",
            wd,
            [
                os.path.join(wd, "reference.fasta"),
                self.state["dbinfo"]["files"]["backbone_tree"],
                os.path.join(wd, "query.fasta"),
            ],
            ["epa_result.jplace"],
            model=self.state["dbinfo"]["iqtree_model"],
            cores=self.state["threads_epa"],
        )
        # Return False if placement failed
        if R.success is False:
            return None
        # load the results
        try:
            with open(os.path.join(wd, "epa_result.jplace"), "r") as info_file:
                data = info_file.read()
                info = json.loads(data)
        except json.decoder.JSONDecodeError:
            jfile = os.path.join(wd, "epa_result.jplace")
            logging.error(
                "Malformed epa-ng result file. Did you ungracefully terminate? Remove the file {} and try again".format(
                    jfile
                )
            )
            exit(203)
        info["tog"] = self.guppy_tree(self.state["workdir"], os.path.join(wd, "epa_result.jplace"))
        logging.debug("Placed {} markers in the reference tree".format(len(info["placements"])))
        return info

    def search_markers(self, use_all=False):
        """
        Search for marker genes needed to place the genome in the reference tree

        :param fasta_faa: Fasta file with proteins
        :param workdir: Folder to create temp direcories in

        :return: list with names of marker genes
        """
        wd = os.path.join(self.state["workdir"], "hmm", "placement", self.state["clade"])
        if use_all:
            # for using the entrie db
            hmm_file = self.state["dbinfo"]["files"]["hmm_db"]
        else:
            # normal use case
            hmm_file = self.state["dbinfo"]["files"]["hmm_placement_search"]
        outfile = "found_markers.tbl"
        hmmsearch(
            "hmmsearch",
            workdir=wd,
            infiles=[hmm_file, self.state["faa"]],
            outfiles=[outfile],
            cores=self.state["threads"],
            cut_ga=True,
        )
        data = read_hmmer(os.path.join(wd, outfile))
        file.delete_but(wd, outfile)
        return data

    def load_db(self, db, clade="base"):
        """
        Function to check if the database has all the files we expect

        :param db: Path to database folder
        """
        logging.debug("Validating db")
        info = self._load_db_info(db, "db_" + clade)
        if "files" not in info.keys():
            logging.error("Malformed info file")
            raise ValueError("The database seems to be malformed")

        def rec_check_path(d):
            for key, p in d.items():
                if type(p) is str:
                    if not os.path.exists(p) and not os.path.exists("{}.h3f".format(p)):
                        raise FileNotFoundError("Database is missing files: {}".format(p))
                elif type(p) is dict:
                    rec_check_path(p)
            return d

        # for each file, make sure we can find it
        rec_check_path(info["files"])
        self.state["loaded_dbs"][clade] = info

        return info

    def _load_db_info(self, db, clade="base"):
        """
        Load required information from the db

        :param db: Path to database folder
        :return: json from database info file
        """
        # read in metadata from json file
        logging.debug("Loading database info for subdb {}".format(clade))
        info_file = os.path.join(db, clade, "db_info.json")
        if not os.path.exists(info_file):
            raise FileNotFoundError("Could not find database info file")

        # read file
        with open(info_file, "r") as info_file:
            data = info_file.read()
            info = json.loads(data)

        def rec_add_path(d, prefix=""):
            for key, val in d.items():
                if type(val) is str:
                    d[key] = os.path.join(prefix, val)
                elif type(val) is dict:
                    rec_add_path(val, prefix)
            return d

        # recursivly add prefix of db location
        rec_add_path(info["files"], os.path.join(db, clade))

        # placement profiles we load from the folder
        place_folder = info["files"]["hmm_placement_place"]
        info["files"]["hmm_placement_place"] = {}
        for f in os.listdir(place_folder):
            n = os.path.basename(f)
            profile = n.rsplit(".", 1)[0]
            suffix = n.rsplit(".", 1)[1]
            if profile not in info["files"]["hmm_placement_place"].keys():
                info["files"]["hmm_placement_place"][profile] = {}
            info["files"]["hmm_placement_place"][profile][suffix] = os.path.realpath(os.path.join(place_folder, f))

        info["dblocation"] = os.path.abspath(db)
        return info

    def check_dependencies(self):
        """
        check that all dependencies for a successfull run
        are accesible.

        :param self.state['seqtype']: either DNA or AA
        """
        logging.debug("Checking dependencies")
        base_dep = ["epa-ng", "hmmsearch", "hmmalign", "hmmfetch", "guppy"]
        pred_dep = ["metaeuk"]

        if self.state["seqtype"] == "DNA":
            base_dep.extend(pred_dep)

        missing = []
        for dep in base_dep:
            logging.debug("Checking dependency {}".format(dep))
            if which(dep) is None:
                missing.append(dep)

        if len(missing) > 0:
            logging.error("Missing dependencies.")
            raise EnvironmentError("Missing dependencies: {}".format(", ".join(missing)))

    def write_result(self, result, outfile="eukcc.tsv", sep="\t"):
        """
        Given a dictionary will write a nice column output file.
        If a list is given will write a row for each entry, assuming each
        element is a dictionary with the right fields
        """
        if result is None:
            logging.debug("Empty result given, so we write NAs for all values")

        if type(result) == dict:
            result = [result]
        elif type(result) is eukcc_state:
            self.state = result
            result = [self.state]
        # file.isdir(os.path.dirname(outfile))
        fields = ["fasta", "completeness", "contamination"]
        with open(outfile, "w") as fout:
            writer = csv.DictWriter(fout, fields, delimiter=sep, extrasaction="ignore")
            writer.writeheader()
            for res in result:
                row = {
                    "fasta": res["fasta"],
                    "completeness": res["quality"]["completeness"],
                    "contamination": res["quality"]["contamination"],
                }
                writer.writerow(row)

        logging.info("Wrote output to: {}".format(outfile))
        self.state.save_state()

    def write_extra(self, state):
        """
        writes some extra tables for the output
        """
        if "scmgs_table" in state.keys():
            with gzip.open(os.path.join(state["out"], "scmg_marker_table.csv.gz"), "wt") as fout:
                c = csv.DictWriter(fout, delimiter="\t", fieldnames=state["scmgs_table"][0].keys())
                c.writeheader()
                for row in state["scmgs_table"]:
                    c.writerow(row)
        return

    def load_result(self, path):
        """
        Loads a json.json.gz result file
        """
        if not os.path.exists(path):
            logging.error("Could not find file {}".format(path))
            exit(200)

        with gzip.open(path) as fin:
            data = json.loads(fin.read().decode("utf-8"))
        return data
