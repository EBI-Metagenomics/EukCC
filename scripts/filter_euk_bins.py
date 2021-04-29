#!/usr/bin/env python3
#
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
# provides all file operation functions
# used inthis package
import os
import argparse
import subprocess
import logging
import tempfile
import gzip
from multiprocessing import Pool


# backup fasta handler, so we can use readonly directories
class fa_class:
    def __init__(self, seq, name, long_name):
        self.seq = seq
        self.name = name
        self.long_name = long_name

    def __str__(self):
        return self.seq

    def __len__(self):
        return len(self.seq)


def Fasta(path):
    """
    Iterator for fasta files
    """
    entry = False

    with open(path) as fin:
        for line in fin:
            if line.startswith(">"):
                if entry is not False:
                    entry.seq = "".join(entry.seq)
                    yield entry
                # define new entry
                long_name = line.strip()[1:]
                name = long_name.split()[0]
                entry = fa_class([], name, long_name)
            else:
                entry.seq.append(line.strip())
        # yield last one
        entry.seq = "".join(entry.seq)
        yield entry


def gunzip(path, tmp_dir):
    """
    Gunzip a file for EukRep
    """
    if path.endswith(".gz"):
        fna_path = os.path.join(tmp_dir, "contigs.fna")
        logging.debug("Going to unzip fasta into {}".format(fna_path))
        with gzip.open(path, "r") as fin, open(fna_path, "w") as fout:
            for line in fin:
                fout.write(line.decode())
        path = fna_path
        logging.debug("Done unzipping {}".format(fna_path))

    return path


class EukRep:
    """Class to call and handle EukRep data"""

    def __init__(self, fasta, eukout, bacout=None, minl=1500, tie="euk"):
        self.fasta = fasta
        self.eukout = eukout
        self.bacout = bacout
        self.minl = minl
        self.tie = tie

    def run(self):
        # command list will be called
        cmd = [
            "EukRep",
            "--min",
            str(self.minl),
            "-i",
            self.fasta,
            "--seq_names",
            "-ff",
            "--tie",
            self.tie,
            "-o",
            self.eukout,
        ]
        if self.bacout is not None:
            cmd.extend(["--prokarya", self.bacout])

        subprocess.run(cmd, check=True, shell=False)

        self.read_result()

    def read_result(self):
        self.euks = self.read_eukfile(self.eukout)

        self.bacs = set()
        if self.bacout is not None:
            self.bacs = self.read_eukfile(self.bacout)

    def read_eukfile(self, path):
        lst = set()
        with open(path) as infile:
            for line in infile:
                lst.add(line.strip())
        return lst


class bin:
    def __init__(self, path, eukrep):
        self.e = eukrep
        self.bname = os.path.basename(path)
        self.path = os.path.abspath(path)

    def __str__(self):
        return "{} euks: {} bacs: {}".format(self.bname, self.table["euks"], self.table["bacs"])

    def stats(self):
        """read bin content and figure genomic composition"""
        logging.debug("Loading bin")
        fa_file = Fasta(self.path)
        stats = {"euks": 0, "bacs": 0, "NA": 0, "sum": 0}
        # loop and compute stats
        logging.debug("Make per bin stats")
        for seq in fa_file:
            if seq.name in self.e.euks:
                stats["euks"] += len(seq)
            elif seq.name in self.e.bacs:
                stats["bacs"] += len(seq)
            else:
                stats["NA"] += len(seq)

        stats["sum"] = sum([v for k, v in stats.items()])

        self.table = stats

    def decide(self, eukratio=0.2, bacratio=0.1, minbp=100000, minbpeuks=1000000):
        """
        rule to handle decision logic
        """
        keep = True
        allb = self.table["sum"]
        if self.table["euks"] < minbpeuks:
            keep = False
            logging.info(f"Eukaryotic DNA amount only {self.table['euks']} instead of target {minbpeuks}")

        elif self.table["euks"] / allb <= eukratio:
            keep = False
            logging.info(f"Eukaryotic DNA ratio not higher than {eukratio}")

        elif self.table["bacs"] / allb >= bacratio:
            keep = False
            logging.info(f"Bacterial DNA content higher than {bacratio}")

        elif self.table["sum"] < minbp:
            keep = False
            logging.info("We did not find at least %d bp of DNA", minbp)

        self.keep = keep


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", help="path for the output table", default="assignment.csv", type=str)
    parser.add_argument("bins", nargs="+", help="all bins to classify", type=str)
    parser.add_argument(
        "--threads",
        "-t",
        type=int,
        help="How many bins should be run in parallel (Default: 1)",
        default=1,
    )
    parser.add_argument(
        "--minl",
        type=int,
        help="define minimal length of contig for EukRep \
                        to classify (default: 1500)",
        default=1500,
    )
    parser.add_argument(
        "--eukratio",
        type=float,
        help="This ratio of eukaryotic DNA to all DNA has to be found\
                    at least (default: 0, ignore)",
        default=0,
    )
    parser.add_argument(
        "--bacratio",
        type=float,
        help="discard bins with bacterial ratio of higher than\
                    (default: 1, ignore)",
        default=1,
    )
    parser.add_argument(
        "--minbp",
        type=float,
        help="Only keep bins with at least n bp of dna\
                    (default: 8000000)",
        default=8000000,
    )
    parser.add_argument(
        "--minbpeuks",
        type=float,
        help="Only keep bins with at least n bp of Eukaryotic dna\
                    (default: 5000000)",
        default=5000000,
    )
    parser.add_argument("--rerun", action="store_true", help="rerun even if output exists", default=False)
    parser.add_argument("--quiet", action="store_true", help="supress information", default=False)
    parser.add_argument("--debug", action="store_true", help="Make it more verbose", default=False)

    args = parser.parse_args()

    # define logging
    logLevel = logging.INFO
    if args.quiet:
        logLevel = logging.WARNING
    elif args.debug:
        logLevel = logging.DEBUG
    logging.basicConfig(format="%(asctime)s %(message)s", datefmt="%m/%d/%Y %H:%M:%S: ", level=logLevel)

    def evaluate_bin(path):
        if not os.path.exists(path):
            logging.error("Can not find file {}".format(path))
            exit(1)

        logging.info("Launch on {}".format(path))
        with tempfile.TemporaryDirectory(prefix="filter_EukRep_") as tmp_dir:
            logging.debug("Using tmp folder: {}".format(tmp_dir))
            eukfile = os.path.join(tmp_dir, "euks.fna")
            bacfile = os.path.join(tmp_dir, "bacs.fna")
            # EukRep can not deal with Gzipped Fasta files, so we unzip it in case it is a Gzip file
            path = gunzip(path, tmp_dir)
            # Launching EukRep
            logging.debug(f"Starting EukRep on {path}")
            eukrep_result = EukRep(path, eukfile, bacfile, minl=args.minl)
            eukrep_result.run()
            b = bin(path, eukrep_result)
            b.stats()
            b.decide(eukratio=args.eukratio, bacratio=args.bacratio, minbp=args.minbp, minbpeuks=args.minbpeuks)
            return b

    # multithreading pool
    pool = Pool(processes=args.threads)
    results = pool.map(evaluate_bin, args.bins)
    pool.close()
    pool.join()

    with open(args.output, "w") as outfile:
        outfile.write("path,binname,passed,bp_eukaryote,bp_prokaryote,bp_unassigned,bp_sum\n")
        for b in results:
            outfile.write(
                f"{b.path},{b.bname},{b.keep},{b.table['euks']},{b.table['bacs']},{b.table['NA']},{b.table['sum']}\n"
            )
