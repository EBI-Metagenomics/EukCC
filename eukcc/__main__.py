import os
import logging
from types import SimpleNamespace
import click
from eukcc import __version__
from eukcc.eukcc import eukcc, eukcc_state
from eukcc.fasta import determine_type
from eukcc.file import file
from eukcc.refine import eukcc_folder
from eukcc.set_define import define_set_w_genomes


CONTEXT_SETTINGS = {"help_option_names": ["-h", "--help"]}


def setup_logging(quiet, debug, out=None):
    """
    Log to the terminal and, if an output folder is given, to <out>/eukcc.log
    """
    logLevel = logging.INFO
    if quiet:
        logLevel = logging.WARNING
    elif debug:
        logLevel = logging.DEBUG

    handlers = [logging.StreamHandler()]
    if out is not None:
        # create output folder if not exists
        file.isdir(out)
        handlers.insert(0, logging.FileHandler(os.path.join(out, "eukcc.log")))

    logging.basicConfig(
        format="%(asctime)s %(message)s",
        datefmt="%d-%m-%Y %H:%M:%S: ",
        level=logLevel,
        handlers=handlers,
    )


@click.group(
    context_settings=CONTEXT_SETTINGS,
    invoke_without_command=True,
    help="Framework to evaluate completeness and contamination of a eukaryotic MAG or isolate genome. Only valid for microbial eukaryotes.",
)
@click.option("--quiet", "-q", is_flag=True, default=False, help="Silence most output")
@click.option("--debug", "-d", is_flag=True, default=False, help="Debug and thus ignore safety")
@click.version_option(__version__, "-v", "--version", message="EukCC version %(version)s")
@click.pass_context
def main(ctx, quiet, debug):
    # without a subcommand, show the help and exit with 0
    if ctx.invoked_subcommand is None:
        click.echo(ctx.get_help())
        ctx.exit(0)
    ctx.obj = {"quiet": quiet, "debug": debug}


# ----------------------------------------------------------------------------
# eukcc single
# ----------------------------------------------------------------------------

publications = {
    "eukcc": 'Saary, Paul, Alex L. Mitchell, and Robert D. Finn. "Estimating the quality of eukaryotic genomes recovered from metagenomic analysis with EukCC." Genome biology 21.1 (2020): 1-21.',
    "hmmer": 'Eddy, Sean R. "Accelerated profile HMM searches." PLoS Comput Biol 7.10 (2011): e1002195.',
    "pplacer": 'Matsen, Frederick A., Robin B. Kodner, and E. Virginia Armbrust. "pplacer: linear time maximum-likelihood and Bayesian phylogenetic placement of sequences onto a fixed reference tree." BMC bioinformatics 11.1 (2010): 1-16.',
    "epa-ng": 'Barbera, Pierre, et al. "EPA-ng: massively parallel evolutionary placement of genetic sequences." Systematic biology 68.2 (2019): 365-369.',
    "metaEuk": 'Levy Karin, Eli, Milot Mirdita, and Johannes Söding. "MetaEuk—Sensitive, high-throughput gene discovery, and annotation for large-scale eukaryotic metagenomics." Microbiome 8 (2020): 1-15.',
}


def no_estimate(E, state):
    """
    It was not possible to estimate this genome. Write the result file with NAs and stop.
    """
    logging.warning("Could not estimate quality, writing NA results")
    E.write_result(E.state, os.path.join(state["out"], "eukcc.tsv"))
    E.terminate()


def run_eukcc(args):
    state = eukcc_state(workdir=None, options=vars(args))
    logging.info(f"EukCC version {__version__}")
    logging.warning(
        "#####################################\nIf you publish using EukCC please make sure to also cite:"
    )
    for software, publication in publications.items():
        logging.warning(f"\n{software}:\n{publication}\n")
    logging.warning("####################################")

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

    if state["max_set_size"] < state["set_size"]:
        logging.error("Max set size cant be smaller than min set size")
        exit(1)

    # validate fasta input
    if state["seqtype"] is None:
        state["seqtype"] = determine_type(state["fasta"])
        logging.info(f"Set sequence type to {state['seqtype']}")

    # launch new EukCC instance
    E = eukcc(state)

    if E.state["seqtype"] == "DNA":
        # predict proteins
        E.predict_protein()
    else:
        # copy fna to faa state entry
        E.state["faa"] = E.state["fasta"]

    if os.stat(E.state["faa"]).st_size == 0:
        logging.warning("No proteins to search for marker genes")
        no_estimate(E, state)

    if E.state["comparison_genomes"] is not None:
        E.state["simple"] = True

    if (
        not E.state["simple"]
        and E.state["clade"] == "base"
        and E.state["use_placement"] is None
    ):
        logging.info("Doing a first pass placement in the base database")
        if E.placement() is None:
            no_estimate(E, state)
        # decide which db to use
        clade = E.determine_subdb()
        if clade != "base":
            E.state["clade"] = clade
            E.state["dbinfo"] = E.load_db(E.state["db"], clade=clade)
        # now that we loaded the new DB we can continue using the rest of the algorythm

    # pick marker set and placement in one step
    if E.pick_marker_set() is None:
        no_estimate(E, state)
    E.state["scmg_data"] = E.hmmsearch_scmg(
        E.state["workdir"], E.state["faa"], E.state["marker_set"]["profiles"]
    )

    logging.debug("Estimating quality")
    E.compute_quality(E.state["scmg_data"], E.state["marker_set"]["profiles"])
    E.add_taxid()

    # aggregate all the interesting data, such as estimated quality
    # the expected and found SCMGs and the workdir
    E.state.save_state()
    E.write_result(E.state, os.path.join(state["out"], "eukcc.tsv"))

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


@main.command(context_settings=CONTEXT_SETTINGS, help="Estimate the quality of a single genome or MAG.")
@click.argument("fasta", type=str)
@click.option("--out", "-o", type=str, default=".", help="Path to output folder (Default: .)")
@click.option("--db", type=str, default=None, help="Path to EukCC DB")
@click.option("--threads", "-t", type=int, default=1, help="Number of threads to use (Default: 1)")
@click.option(
    "--threads_epa",
    type=int,
    default=1,
    help="Number of threads to use for epa-ng, recommended: 1 (Default: 1)",
)
@click.option("--set_size", type=int, default=20, help="Minimal number of marker genes to use (Default: 20)")
@click.option(
    "--max_set_size",
    type=int,
    default=500,
    help="Maximal number of marker genes used, set to 0 to include all possible marker genes (Default: 500)",
)
@click.option(
    "--marker_prevalence",
    type=float,
    default=95,
    help="Percentage of species in which markers should be found (Default: 95)",
)
@click.option("--DNA", "seqtype", flag_value="DNA", default=None, help="The fasta file contains DNA sequences")
@click.option("--AA", "seqtype", flag_value="AA", help="The fasta file contains amino acid sequences")
@click.option(
    "--taxids",
    type=str,
    multiple=True,
    help="Taxid to use as set starting point, repeat the option for several taxids",
)
@click.option(
    "--genomes",
    "comparison_genomes",
    type=str,
    multiple=True,
    help="Genome file to base a SCMG set upon, repeat the option for several genomes",
)
@click.option(
    "--use_placement",
    type=str,
    default=None,
    help="Path to previous result file, to use exact same marker gene set",
)
@click.option(
    "--set_number_species",
    type=int,
    default=3,
    help="Minimal number of species to define a set. Reduce this if no sets can be found (Default: 3)",
)
@click.option(
    "--select_best_guess",
    "set_selection",
    flag_value="best_guess",
    default="best_guess",
    help="Use best guess to select marker gene set (Default)",
)
@click.option(
    "--select_species",
    "set_selection",
    flag_value="species",
    help="Use species count to select best marker gene set (Default: best guess)",
)
@click.option(
    "--use_ncbi_tree",
    "use_ncbi",
    is_flag=True,
    default=False,
    help="Instead of using the EukCC phylogenetic tree, rely on NCBI taxids (default: False)",
)
@click.option(
    "--ignore_tree",
    is_flag=True,
    default=False,
    help="Advanced option, mainly for debugging. Can ignore the tree if genomes are known via taxids for example",
)
@click.option(
    "--simple",
    is_flag=True,
    default=False,
    help="Use global DB instead of clade specific dbs (faster, not suitable for protozoa)",
)
@click.option("--clade", type=str, default="base", help="Define clade as base, fungi, protozoa or plant")
@click.option(
    "--rerun",
    "-r",
    is_flag=True,
    default=False,
    help="Rerun and remove any previously computed data in the target folder",
)
@click.option(
    "--no_dynamic_root",
    is_flag=True,
    default=False,
    help="Do not re-root tree dynamically, to get best set detection (default: True)",
)
@click.option("--extra", is_flag=True, default=False, help="Produce extra output files.")
@click.option("--keep", is_flag=True, default=False, help="Keep workdir after the run.")
@click.pass_obj
def single(obj, **kwargs):
    setup_logging(obj["quiet"], obj["debug"], kwargs["out"])
    kwargs["taxids"] = list(kwargs["taxids"]) or None
    kwargs["comparison_genomes"] = list(kwargs["comparison_genomes"]) or None
    kwargs["dynamic_root"] = not kwargs.pop("no_dynamic_root")
    kwargs["clean"] = not kwargs.pop("keep")
    run_eukcc(SimpleNamespace(**obj, **kwargs))


# ----------------------------------------------------------------------------
# eukcc folder
# ----------------------------------------------------------------------------

@main.command(context_settings=CONTEXT_SETTINGS, help="Estimate the quality of a folder of bins and refine them.")
@click.argument("binfolder", type=str)
@click.option("--out", "-o", type=str, default=".", help="Path to output folder (Default: .)")
@click.option("--db", type=str, default=None, help="Path to EukCC DB")
@click.option("--threads", "-t", type=int, default=1, help="Number of threads to use (Default: 1)")
@click.option(
    "--threads_epa",
    type=int,
    default=1,
    help="Number of threads to use for epa-ng, recommended: 1 (Default: 1)",
)
@click.option(
    "--marker_prevalence",
    type=float,
    default=95,
    help="Percentage of species in which markers should be found (Default: 95)",
)
@click.option(
    "--links",
    type=str,
    default=None,
    help="Path to a link table generated with binlinks.py. If supplied, paired reads will be used to refine bins (Recommended)",
)
@click.option(
    "--min_links",
    type=int,
    default=100,
    help="Number of paired reads matching between bins for merging to happen (default: 100)",
)
@click.option("--prefix", type=str, default="merged.", help="Prefix to add for merged bins (default: merged.)")
@click.option(
    "--improve_ratio",
    type=float,
    default=5,
    help="Ratio of completeness to contamination change (Default: 5)",
)
@click.option(
    "--improve_percent",
    type=float,
    default=10,
    help="A merger must increase completeness at least by n percent (Default: 10)",
)
@click.option(
    "--n_combine",
    type=int,
    default=1,
    help="How many small bins should be merged into a medium sized bin (Default: 1)",
)
@click.pass_obj
def folder(obj, **kwargs):
    setup_logging(obj["quiet"], obj["debug"], kwargs["out"])
    eukcc_folder(SimpleNamespace(**obj, **kwargs))


# ----------------------------------------------------------------------------
# eukcc define_set
# ----------------------------------------------------------------------------

@main.command("define_set", context_settings=CONTEXT_SETTINGS, help="Define a marker gene set from a group of genomes.")
@click.argument("genomes", type=str, nargs=-1, required=True)
@click.option("--out", "-o", type=str, default=".", help="Path to output folder (Default: .)")
@click.option("--db", type=str, default=None, help="Path to EukCC DB")
@click.option("--threads", "-t", type=int, default=1, help="Number of threads to use (Default: 1)")
@click.option("--set_size", type=int, default=20, help="Minimal number of marker genes to use (Default: 20)")
@click.option(
    "--max_set_size",
    type=int,
    default=500,
    help="Maximal number of marker genes used, set to 0 to include all possible marker genes (Default: 500)",
)
@click.option(
    "--marker_prevalence",
    type=float,
    default=95,
    help="Percentage of species in which markers should be found (Default: 95)",
)
@click.option("--name", "-n", type=str, default="set", help="Name of the set (default: set)")
@click.option("--DNA", "seqtype", flag_value="DNA", default=None, help="The fasta file contains DNA sequences")
@click.option(
    "--AA",
    "seqtype",
    flag_value="AA",
    help="All provided files are proteins, no gene prediction will be performed. Omit if input is mixed or genomic",
)
@click.pass_obj
def define_set(obj, **kwargs):
    setup_logging(obj["quiet"], obj["debug"], kwargs["out"])
    kwargs["genomes"] = list(kwargs["genomes"])
    define_set_w_genomes(SimpleNamespace(**obj, **kwargs))


# ----------------------------------------------------------------------------
# eukcc ncbi_update
# ----------------------------------------------------------------------------

@main.command("ncbi_update", context_settings=CONTEXT_SETTINGS, help="Update the NCBI taxonomy database used by ete3.")
@click.pass_obj
def ncbi_update(obj):
    from ete3 import NCBITaxa

    setup_logging(obj["quiet"], obj["debug"])
    logging.info("Going to fetch NCBI info using ete3")
    NCBITaxa().update_taxonomy_database()
