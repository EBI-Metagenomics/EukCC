# This file is part of the EukCC (https://github.com/openpaul/eukcc).
# Copyright (c) 2019 Paul Saary
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

import logging
import configargparse
from eukcc import workflow
import os


def main():
    # set arguments
    # arguments are passed to classes
    parser = configargparse.ArgumentParser(description="Evaluate completeness and contamination of a MAG.")
    parser.add_argument("fasta", type=str, help="Run script on this bin (fasta file)")
    parser.add_argument("--db", type=str, required=True, help="Path to EukCC DB")
    parser.add_argument(
        "--outdir",
        "-o",
        type=str,
        default="./",
        help="Location for the output. Names will be prefixed using the bin filenames",
    )
    parser.add_argument(
        "--config", "-c", type=str, required=False, is_config_file=True, help="Config file to define parameters, YAML",
    )
    parser.add_argument(
        "--ncores",
        "-n",
        metavar="int",
        type=int,
        default=1,
        help="set number of cores for GeneMark-ES, pplacer and Hmmer",
    )
    parser.add_argument(
        "--ncorespplacer",
        metavar="int",
        type=int,
        default=0,
        help="Pplacer requires a lot of memory. If you want \
                              you can set less cores for pplacer,\
                              which improves memory consumption significantly",
    )
    parser.add_argument(
        "--hmm", dest="hmm", type=str, default=None, help="run hmmer on all these HMMs instead",
    )
    parser.add_argument(
        "--training",
        dest="training",
        action="store_true",
        default=False,
        help="Run EukCC in training mode (needed to create a new release of the DB)",
    )
    parser.add_argument(
        "--bed",
        "-b",
        metavar="file.bed",
        type=str,
        default=None,
        help="Pass bedfile if you called genes manually. \
                        Assumes only a single fasta (protein) is passed and implies --noglob",
    )
    parser.add_argument(
        "--force",
        "-f",
        dest="force",
        action="store_true",
        default=False,
        help="Force rerun of computation even if \
                                              output is newer than input. Don't resume previous run.",
    )
    parser.add_argument(
        "--keeptemp",
        dest="clean",
        action="store_false",
        default=True,
        help="Keep all temporary files, by default EukCC will remove some temp files",
    )
    parser.add_argument(
        "--fplace",
        "-p",
        dest="fplace",
        action="store_true",
        default=False,
        help="Force rerun of placement and subsequent steps",
    )
    parser.add_argument(
        "--noglob", "-g", dest="noglob", action="store_true", default=False, help="Do not expand paths using glob",
    )
    parser.add_argument(
        "--quiet", "-q", dest="quiet", action="store_true", default=False, help="Silcence most output",
    )
    parser.add_argument(
        "--debug", "-d", action="store_true", default=False, help="Debug and thus ignore safety",
    )
    parser.add_argument(
        "--HPA", default=False, action="store_true", help="Set placement method to HPA",
    )
    parser.add_argument(
        "--nPlacements",
        type=int,
        default=2,
        metavar="n",
        help="Set number of proteins to support location \
                                in tree (default: 2)",
    )
    parser.add_argument(
        "--minGenomes",
        type=int,
        default=3,
        metavar="n",
        help="Minimal number of genomes to support a set (default: 3)",
    )
    parser.add_argument(
        "--fullineage", default=False, action="store_true", help="Output full lineage for MAGs",
    )
    parser.add_argument(
        "--minPlacementLikelyhood",
        default=0.4,
        type=float,
        metavar="float",
        help="minimal pplacer likelyhood (default: 0.4)",
    )
    parser.add_argument(
        "--mindist", type=int, default=2000, metavar="n", help="Distance to collapse hits (default: 2000)",
    )
    parser.add_argument(
        "--touch", default=False, action="store_true", help="Do not run, but touch all output files",
    )
    parser.add_argument(
        "--gmes", default=False, action="store_true", help="only run GeneMark-ES",
    )
    parser.add_argument("--plot", default=False, action="store_true", help="produce plots")
    options = parser.parse_args()

    # define logging
    logLevel = logging.INFO
    if options.quiet:
        logLevel = logging.WARNING
    elif options.debug:
        logLevel = logging.DEBUG
    logging.basicConfig(
        format="%(asctime)s %(message)s", datefmt="%m/%d/%Y %H:%M:%S: ", level=logLevel,
    )
    logging.debug("Launching EukCC in debug mode")
    logging.info("Starting EukCC")

    # Now we start the run with EukCC
    # All magic numbers should be defined in info.py if they are not
    # part of the configuration options
    m = workflow.eukcc(options)

    # skip gene predition if this is already protein sequences
    if options.bed is None:
        # run gmes
        proteinfaa, bedfile = m.gmes(options.fasta)
    else:
        proteinfaa = options.fasta
        bedfile = options.bed

    # terminate if only gmes step was to be run
    if m.cfg["gmes"]:
        logging.info("Finished running GeneMark-ES")
        logging.info("Terminating as requested")
        exit(0)

    # run hmm file if we are asked to
    # this is needed during for training
    if m.cfg["training"] or m.cfg["hmm"]:
        logging.info("Running on custom hmm for training mode")
        m.runPlacedHMM(m.cfg["hmm"], proteinfaa, bedfile)
        logging.info("Stopping now as we are only doing training")
        exit()

    # place using pplacer and hmmer
    m.place(proteinfaa, bedfile)

    # concat hmms for hmmer
    hmmfile = m.concatHMM()
    # run Hmmer for sets of placement
    hits = m.runPlacedHMM(hmmfile, proteinfaa, bedfile)
    # infer lineage
    _ = m.inferLineage(m.placements[m.cfg["placementMethod"]])

    # estimate completeness and contamiantion
    outputfile = os.path.join(m.cfg["outdir"], "eukcc.tsv")
    m.estimate(hits, outputfile, m.placements[m.cfg["placementMethod"]])

    if m.cfg["plot"]:
        _ = m.plot()
