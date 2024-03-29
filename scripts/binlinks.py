#!/usr/bin/env python3
import pysam
from Bio import SeqIO
from collections import defaultdict
import os
import argparse
import logging
import csv


def is_in(read, contig_map, within=1000):
    if read.reference_name not in contig_map.keys():
        return False
    if read.reference_start <= within or read.reference_end <= within:
        return True
    elif read.reference_start > (
        contig_map[read.reference_name] - within
    ) or read.reference_end > (contig_map[read.reference_name] - within):
        return True
    else:
        return False


def keep_read(read, contig_map, within=1000, min_ANI=98, min_cov=0):
    ani = (
        (read.query_alignment_length - read.get_tag("NM"))
        / float(read.query_alignment_length)
        * 100
    )
    cov = read.query_alignment_length / float(read.query_length) * 100

    if ani >= min_ANI and cov >= min_cov and is_in(read, contig_map, within) is True:
        return True
    else:
        return False


def contig_map(bindir, suffix=".fa"):
    m = {}
    for f in os.listdir(bindir):
        if f.endswith(suffix) is False:
            continue
        path = os.path.join(bindir, f)
        with open(path, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                m[record.name] = len(record.seq)
    return m


def bin_map(bindir, suffix=".fa"):
    contigs = defaultdict(str)
    contigs_per_bin = defaultdict(int)
    for f in os.listdir(bindir):
        if f.endswith(suffix) is False:
            continue
        path = os.path.join(bindir, f)
        binname = os.path.basename(f)
        with open(path, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                contigs[record.name] = binname
                contigs_per_bin[binname] += 1
    return contigs, contigs_per_bin


def read_pair_generator(bam):
    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    From: https://www.biostars.org/p/306041/
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch():
        if not read.is_paired or read.is_secondary or read.is_supplementary:
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]


def read_bam_file(bamf, link_table, cm, within, ANI):
    samfile = pysam.AlignmentFile(bamf, "rb")

    # generate link table
    logging.info("Parsing Bam file. This can take a few moments")
    for read, mate in read_pair_generator(samfile):
        if keep_read(read, cm, within, min_ANI=ANI) and keep_read(
            mate, cm, within, min_ANI=ANI
        ):
            # fill in the table
            link_table[read.reference_name][mate.reference_name] += 1
            if read.reference_name != mate.reference_name:
                link_table[mate.reference_name][read.reference_name] += 1

    return link_table


def main():
    # set arguments
    # arguments are passed to classes
    parser = argparse.ArgumentParser(
        description="Evaluate completeness and contamination of a MAG."
    )
    parser.add_argument("bindir", type=str, help="Run script on these bins")
    parser.add_argument(
        "bam",
        type=str,
        help="Bam file(s) with reads aligned against all contigs making up the bins",
        nargs="+",
    )
    parser.add_argument(
        "--out",
        "-o",
        type=str,
        required=False,
        help="Path to output table (Default: links.csv)",
        default="links.csv",
    )
    parser.add_argument(
        "--ANI", type=float, required=False, help="ANI of matching read", default=99
    )
    parser.add_argument(
        "--within",
        type=int,
        required=False,
        help="Within this many bp we need the read to map",
        default=1000,
    )
    parser.add_argument(
        "--contigs",
        "-c",
        action="store_true",
        default=False,
        help="Instead of bins print contigs",
    )
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
    args = parser.parse_args()

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
    )

    bindir = args.bindir

    cm = contig_map(bindir)
    bm, contigs_per_bin = bin_map(bindir)
    logging.debug("Found {} contigs".format(len(cm)))

    link_table = defaultdict(lambda: defaultdict(int))
    bin_table = defaultdict(lambda: defaultdict(int))

    # iterate all bam files
    for bamf in args.bam:
        link_table = read_bam_file(bamf, link_table, cm, args.within, args.ANI)

    logging.debug("Created link table with {} entries".format(len(link_table)))

    # generate bin table
    for contig_1, dic in link_table.items():
        for contig_2, links in dic.items():
            bin_table[bm[contig_1]][bm[contig_2]] += links

    logging.debug("Created bin table with {} entries".format(len(bin_table)))

    out_data = []
    logging.debug("Constructing output dict")
    if args.contigs:
        for contig_1, linked in link_table.items():
            for contig_2, links in linked.items():
                out_data.append(
                    {
                        "bin_1": bm[contig_1],
                        "bin_2": bm[contig_2],
                        "contig_1": contig_1,
                        "contig_2": contig_2,
                        "links": links,
                        "bin_1_contigs": contigs_per_bin[bm[contig_1]],
                        "bin_2_contigs": contigs_per_bin[bm[contig_2]],
                    }
                )
    else:
        for bin_1, dic in bin_table.items():
            for bin_2, links in dic.items():
                out_data.append({"bin_1": bin_1, "bin_2": bin_2, "links": links})

    logging.debug("Out data has {} rows".format(len(out_data)))
    # results
    logging.info("Writing output")
    with open(args.out, "w") as fout:
        if len(out_data) > 0:
            cout = csv.DictWriter(fout, fieldnames=list(out_data[0].keys()))
            cout.writeheader()
            for row in out_data:
                cout.writerow(row)
        else:
            logging.warning("No rows to write")


if __name__ == "__main__":
    main()
