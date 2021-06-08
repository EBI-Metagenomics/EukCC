import os
import argparse
import logging
from eukcc.eukcc import eukcc, eukcc_state
from eukcc.fasta import determine_type
import eukcc.version as version
from eukcc.refine import eukcc_folder


publications = {
    "eukcc": 'Saary, Paul, Alex L. Mitchell, and Robert D. Finn. "Estimating the quality of eukaryotic genomes recovered from metagenomic analysis with EukCC." Genome biology 21.1 (2020): 1-21.',
    "hmmer": 'Eddy, Sean R. "Accelerated profile HMM searches." PLoS Comput Biol 7.10 (2011): e1002195.',
    "pplacer": 'Matsen, Frederick A., Robin B. Kodner, and E. Virginia Armbrust. "pplacer: linear time maximum-likelihood and Bayesian phylogenetic placement of sequences onto a fixed reference tree." BMC bioinformatics 11.1 (2010): 1-16.',
    "epa-ng": 'Barbera, Pierre, et al. "EPA-ng: massively parallel evolutionary placement of genetic sequences." Systematic biology 68.2 (2019): 365-369.',
    "metaEuk": 'Levy Karin, Eli, Milot Mirdita, and Johannes Söding. "MetaEuk—Sensitive, high-throughput gene discovery, and annotation for large-scale eukaryotic metagenomics." Microbiome 8 (2020): 1-15.',
}


def update():
    from ete3 import NCBITaxa

    logging.info("Going to fetch NCBI info using ete3")
    ncbi = NCBITaxa()

    ncbi.update_taxonomy_database()


def run_eukcc(args):
    state = eukcc_state(workdir=None, options=vars(args))
    logging.info("EukCC version {}".format(version.__version__))
    logging.warning(
        "#####################################\nIf you publish using EukCC please make sure to also cite:"
    )
    for software, publication in publications.items():
        logging.warning("\n{}:\n{}\n".format(software, publication))
    logging.warning("####################################")

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

    # validate fasta input
    if state["seqtype"] is None:
        state["seqtype"] = determine_type(state["fasta"])
        logging.info("Set sequence type to {}".format(state["seqtype"]))

    # launch new EukCC instance
    E = eukcc(state)

    if E.state["seqtype"] == "DNA":
        # predict proteins
        E.predict_protein()
    else:
        # copy fna to faa state entry
        E.state["faa"] = E.state["fasta"]

    if not E.state["simple"] and E.state["clade"] == "base":
        logging.info("Doing a first pass placement in the base database")
        if E.placement() is None:
            E.terminate(1)
        # decide which db to use
        clade = E.determine_subdb()
        if clade != "base":
            E.state["clade"] = clade
            E.state["dbinfo"] = E.load_db(E.state["db"], clade=clade)
        # now that we loaded the new DB we can continue using the rest of the algorythm

    # pick marker set and placement in one step
    if E.pick_marker_set() is None:
        E.terminate(1)
    E.state["scmg_data"] = E.hmmsearch_scmg(
        E.state["workdir"], E.state["faa"], E.state["marker_set"]["profiles"]
    )

    logging.debug("Estimating quality")
    E.compute_quality(E.state["scmg_data"], E.state["marker_set"]["profiles"])

    # aggregate all the interesting data, such as estimated quality
    # the expected and found SCMGs and the workdir
    E.state.save_state()
    E.write_result(E.state, os.path.join(state["out"], "eukcc.csv"))
    if state["extra"]:
        # search for missing markers

        E.state["scmgs_table"] = E.hmmsearch_scmg(
            state["workdir"],
            state["faa"],
            state["marker_set"]["profiles"],
            cut_ga=False,
        )

        E.write_extra(E.state)

    if state["clean"]:
        E.remove_workdir(state["workdir"])
    E.state.checkpoint("Done")
    return


def main():
    # set arguments
    # arguments are passed to classes
    parser = argparse.ArgumentParser(
        description="Framework to evaluate completeness and contamination of a eukaryotic MAG or isolate genome. Only valid for microbial eukaryotes."
    )

    # define subcommands
    subparsers = parser.add_subparsers(title="Subcommands", dest="command")
    pars_eukcc = subparsers.add_parser("single")
    pars_folder = subparsers.add_parser("folder")
    subparsers.add_parser("ncbi_update")

    parser.add_argument(
        "--quiet",
        "-q",
        dest="quiet",
        action="store_true",
        default=False,
        help="Silcence most output",
    )
    parser.add_argument(
        "--debug",
        "-d",
        action="store_true",
        default=False,
        help="Debug and thus ignore safety",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="EukCC version {}".format(version.__version__),
    )

    pars_folder.add_argument(
        "binfolder", type=str, help="Run script on bins in this folder"
    )
    pars_folder.add_argument(
        "--links",
        type=str,
        required=False,
        help="Path to a link table generated with bamlinks.py. If suuplied paired reads will be used to refine bins (Recommended)",
        default=None,
    )

    pars_folder.add_argument(
        "--min_links",
        type=int,
        required=False,
        help="Number of paired reads matching between bins for merging to happen (default: 100)",
        default=100,
    )
    pars_folder.add_argument(
        "--prefix",
        type=str,
        required=False,
        help="Prefix to add for merged bins (default: merged.)",
        default="merged.",
    )

    pars_folder.add_argument(
        "--suffix",
        type=str,
        required=False,
        help="Suffix (default: .fa)",
        default=".fa",
    )
    pars_folder.add_argument(
        "--out",
        "-o",
        type=str,
        required=False,
        help="Path to output folder (Default: .)",
        default=".",
    )
    pars_folder.add_argument("--db", type=str, default=None, help="Path to EukCC DB")
    pars_folder.add_argument(
        "--improve_ratio",
        type=float,
        help="Ratio of completeness to contamination change (Default: 5)",
        default=5,
    )
    pars_folder.add_argument(
        "--improve_percent",
        type=float,
        help="A merger must increase completeness at least by n percent (Default: 10)",
        default=10,
    )
    pars_folder.add_argument(
        "--n_combine",
        type=int,
        help="How many small bins should be merged into a medium sized bin (Default: 1)",
        default=1,
    )
    pars_folder.add_argument(
        "--threads",
        "-t",
        type=int,
        help="Number of threads to use (Default: 1)",
        default=1,
    )
    pars_folder.add_argument(
        "--threads_epa",
        type=int,
        help="Number of threads to use for epa-ng, recommended: 1 (Default: 1)",
        default=1,
    )
    pars_folder.add_argument(
        "--marker_prevalence",
        type=float,
        required=False,
        help="Percentage of species in which markers should be found (Default: 95)",
        default=95,
    )

    # single fasta
    pars_eukcc.add_argument(
        "fasta", type=str, help="Estimate quality of this single bin (fasta file)"
    )
    pars_eukcc.add_argument(
        "--out",
        "-o",
        type=str,
        required=False,
        help="Path to output folder (Default: .)",
        default=".",
    )
    pars_eukcc.add_argument("--db", type=str, default=None, help="Path to EukCC DB")
    pars_eukcc.add_argument(
        "--threads",
        "-t",
        type=int,
        help="Number of threads to use (Default: 1)",
        default=1,
    )
    pars_eukcc.add_argument(
        "--threads_epa",
        type=int,
        help="Number of threads to use for epa-ng, recommended: 1 (Default: 1)",
        default=1,
    )
    pars_eukcc.add_argument(
        "--DNA",
        dest="seqtype",
        action="store_const",
        const="DNA",
        default=None,
        help="The fasta file contains DNA sequenes",
    )
    pars_eukcc.add_argument(
        "--AA",
        dest="seqtype",
        action="store_const",
        const="AA",
        help="The fasta file contains amino acid sequences",
    )
    pars_eukcc.add_argument(
        "--taxids",
        type=str,
        required=False,
        help="Taxids to use as set starting point",
        default=None,
        nargs="+",
    )
    pars_eukcc.add_argument(
        "--set_size",
        type=int,
        required=False,
        help="Minimal number of marker genes to use (Default: 20)",
        default=20,
    )
    pars_eukcc.add_argument(
        "--use_placement",
        type=str,
        required=False,
        help="Path to previous result file, to use exact same marker gene set",
        default=None,
    )
    pars_eukcc.add_argument(
        "--set_number_species",
        type=int,
        required=False,
        help="Minimal number of species to define a set. Reduce this if no sets can be found (Default: 3)",
        default=3,
    )
    pars_eukcc.add_argument(
        "--marker_prevalence",
        type=float,
        required=False,
        help="Percentage of species in which markers should be found (Default: 95)",
        default=95,
    )
    pars_eukcc.add_argument(
        "--max_set_size",
        type=int,
        required=False,
        help="Maximal number of marker genes used, set to 0 to include all possible marker genes (Default: 500)",
        default=500,
    )
    pars_eukcc.add_argument(
        "--select_best_guess",
        dest="set_selection",
        action="store_const",
        const="best_guess",
        default="best_guess",
        help="Use best guess to select marker gene set (Default)",
    )
    pars_eukcc.add_argument(
        "--select_species",
        dest="set_selection",
        action="store_const",
        const="species",
        default="best_guess",
        help="Use species count to select best marker gene set (Default: best guess)",
    )
    pars_eukcc.add_argument(
        "--use_ncbi_tree",
        dest="use_ncbi",
        action="store_true",
        default=False,
        help="Instead of using the EukCC phylogenetic tree, rely on NCBI taxids (default: False)",
    )
    pars_eukcc.add_argument(
        "--gmes",
        dest="use_gmes",
        action="store_true",
        default=False,
        help="Use GeneMark-ES instead of metaeuk (much slower) (default: False)",
    )
    pars_eukcc.add_argument(
        "--ignore_tree",
        action="store_true",
        default=False,
        help="Advanced option, mainly for debugging. Can ignore the tree if genomes are knwon via taxids for example",
    )
    pars_eukcc.add_argument(
        "--simple",
        dest="simple",
        action="store_true",
        default=False,
        help="Use global DB instead of clade specific dbs (faster, not suitable for protozoa)",
    )
    pars_eukcc.add_argument(
        "--clade",
        default="base",
        type=str,
        help="Define clade as base, fungi, protozoa or plants",
    )
    pars_eukcc.add_argument(
        "--rerun",
        "-r",
        dest="rerun",
        action="store_true",
        default=False,
        help="Rerun and remove any previously computed data in the target folder",
    )
    pars_eukcc.add_argument(
        "--no_dynamic_root",
        dest="dynamic_root",
        action="store_false",
        default=True,
        help="Do not re-root tree dynamically, to get best set detection (default: True)",
    )
    pars_eukcc.add_argument(
        "--extra",
        action="store_true",
        dest="extra",
        default=False,
        help="Produce extra output files.",
    )
    pars_eukcc.add_argument(
        "--keep",
        action="store_false",
        dest="clean",
        default=True,
        help="Keep workdir after the run.",
    )
    args = parser.parse_args()
    if args.command is None:
        parser.print_help()
        exit(0)

    # define logging
    logLevel = logging.INFO
    if args.quiet:
        logLevel = logging.WARNING
    elif args.debug:
        logLevel = logging.DEBUG
    logging.basicConfig(
        format="%(asctime)s %(message)s",
        datefmt="%d-%m-%Y %H:%M:%S: ",
        level=logLevel,
        handlers=[
            logging.FileHandler(os.path.join(args.out, "eukcc.log")),
            logging.StreamHandler(),
        ],
    )

    if args.command == "single":
        run_eukcc(args)
    elif args.command == "folder":
        eukcc_folder(args)
    elif args.command == "ncbi_update":
        update()
    else:
        parser.print_help()
