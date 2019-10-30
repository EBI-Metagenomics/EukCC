#!/usr/bin/python3
import os
from pyfaidx import Fasta
import argparse
import subprocess
import logging


def concatenate_bins(bins, outf):
    """simple function to merge all
    bins into scaffolds file"""
    with open(outf, "w") as outfile:
        for path in bins:
            logging.info(f"Concatenating bin {path}")
            with open(path) as infile:
                for line in infile:
                    outfile.write(line)


def create_dir(d):
    if not os.path.isdir(d):
        try:
            os.makedirs(d)
        except OSError as e:
            logging.warning(f"Could not create dir: {d}\n{e}")


class EukRep():
    """Class to call and handle EukRep data"""
    def __init__(self, fasta, eukout, bacout=None, minl=1500,
                 tie='euk'):
        self.fasta = fasta
        self.eukout = eukout
        self.bacout = bacout
        self.minl = minl
        self.tie = tie

    def run(self):
        # command list will be called
        cmd = ['EukRep', '--min', str(self.minl),
               '-i', self.fasta,
               '--seq_names',
               "-ff",
               "--tie", self.tie,
               '-o', self.eukout]
        if self.bacout is not None:
            cmd.extend(['--prokarya', self.bacout])

        subprocess.run(cmd,  check=True, shell=False)

        self.read_result()

    def read_result(self):
        self.euks = self.read_eukfile(self.eukout)

        self.bacs = []
        if self.bacout is not None:
            self.bacs = self.read_eukfile(self.bacout)

    def read_eukfile(self, path):
        lst = []
        with open(path) as infile:
            for line in infile:
                lst.append(line.strip())
        return(lst)


class bin():
    def __init__(self, path, eukrep):
        self.path = path
        self.e = eukrep

    def stats(self):
        """read bin content and figure genomic composition"""
        logging.debug("Loading bin")
        fa_file = Fasta(self.path)
        stats = {'euks': 0,
                 'bacs': 0,
                 'NA': 0,
                 'sum': 0}
        # loop and compute stats
        logging.debug(f"Make per bin stats ({len(fa_file.keys())} contigs)")
        for seq in fa_file:
            if seq.name in self.e.euks:
                stats['euks'] += len(seq)
            elif seq.name in self.e.bacs:
                stats['bacs'] += len(seq)
            else:
                stats['NA'] += len(seq)

        stats['sum'] = sum([v for k, v in stats.items()])

        self.table = stats

    def decide(self, eukratio=0.2, bacratio=0.1):
        """
        rule to handle decision logic
        """
        keep = True
        allb = self.table['sum']
        if self.table['euks']/allb <= eukratio:
            keep = False
            logging.info(f"Rejecting because eukaryotic DNA ratio not higher than {eukratio}")

        if self.table['bacs']/allb >= bacratio:
            keep = False
            logging.info(f"Rejecting because bacterial DNA content higher than {bacratio}")

        self.keep = keep

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output",
                        help="path for the output table",
                        default="asignment.csv",
                        type=str)
    parser.add_argument("bins", nargs="+",
                        help="all bins to classify",
                        type=str)
    parser.add_argument("--tempdir", type=str,
                        help="Will save temp files of the analysis here",
                        default="tmp")
    parser.add_argument("--minl", type=int,
                        help="define minimal length of contig for EukRep \
                        to classify (default: 1500)",
                        default=1500)
    parser.add_argument("--eukratio", type=float,
            help="This ratio of eukaryotic DNA to all DNA has to be found\
                    at least (default: 0.2)",
                        default=0.2)
    parser.add_argument("--bacratio", type=float,
            help="discard bins with bactrial ratio of higher than\
                    (default: 0.1)",
                        default=0.1)
    parser.add_argument("--rerun", action="store_true",
                        help="rerun even if output exists",
                        default=False)
    parser.add_argument("--quiet", action="store_true",
                        help="supress information",
                        default=False)
    parser.add_argument("--debug", action="store_true",
                        help="Make it more verbose",
                        default=False)

    args = parser.parse_args()

    # define logging
    logLevel = logging.INFO
    if args.quiet:
        logLevel = logging.WARNING
    elif args.debug:
        logLevel = logging.DEBUG
    logging.basicConfig(format='%(asctime)s %(message)s',
                        datefmt='%m/%d/%Y %H:%M:%S: ',
                        level=logLevel)

    # creating tmp dir if not exists
    create_dir(args.tempdir)
    # concatenate fastas as we need single file for EukRep
    contigs = os.path.join(args.tempdir, "contigs.fna")
    logging.info("Concatenating bins")
    if args.rerun or not os.path.exists(contigs):
        concatenate_bins(args.bins, contigs)
    else:
        logging.debug("Reusing concatenated contigs")

    # launch EukRep
    eukfile = os.path.join(args.tempdir, "euks.fna")
    bacfile = os.path.join(args.tempdir, "bacs.fna")
    eukrep_result = EukRep(contigs,
                           eukfile,
                           bacfile,
                           minl=args.minl)

    if args.rerun or not os.path.exists(eukfile):
        logging.info("Running EukRep on concatenated contigs")
        eukrep_result.run()
    else:
        logging.debug("Reusing EukRep run from before")
        eukrep_result.read_result()

    with open(args.output, "w") as outfile:
        outfile.write("path,binname,eukaryotic,eukbases,bacbases,unassigned,sum\n")
        for path in args.bins:
            logging.info(f"Deciding on bin: {path}")
            b = bin(path, eukrep_result)
            b.stats()
            b.decide(eukratio=args.eukratio,
                     bacratio=args.bacratio)
            bname = os.path.basename(path)
            outfile.write(f"{path},{bname},{b.keep},{b.table['euks']},{b.table['bacs']},{b.table['NA']},{b.table['sum']}\n")

    
