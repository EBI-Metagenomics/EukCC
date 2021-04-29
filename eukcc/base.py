import os
import gzip
import csv
import logging
from random import seed, sample
from collections import defaultdict
from ete3 import NCBITaxa
from eukcc.fasta import Fasta


def compute_silent_contamination(fna, faa, hits):
    prot_contigs = set()
    for seq in Fasta(faa):
        contig = seq.name.split("_", 1)[1].rsplit("_", 1)[0]
        prot_contigs.add(contig)

    hit_contigs = set()
    for row in hits:
        contig = row["target"].split("_", 1)[1].rsplit("_", 1)[0]
        hit_contigs.add(contig)

    contigs = {}
    total_bp = 0
    for rec in Fasta(fna):
        contigs[rec.name] = len(rec.seq)
        total_bp += len(rec.seq)

    # count all the stats
    stats = {
        "contigs": len(contigs),
        "contigs_w_hits": 0,
        "bp_w_hits": 0,
        "contigs_w_proteins": 0,
        "total_bp": total_bp,
    }
    for contig, size in contigs.items():
        if contig in prot_contigs:
            stats["contigs_w_proteins"] += 1
        if contig in hit_contigs:
            stats["contigs_w_hits"] += 1
            stats["bp_w_hits"] += size
    return stats


def load_tax_info(path, sep=","):
    """
    Load taxonomic information from two column csv
    With first column giving the taxid and second column
    giving the name of the tree leaf.


    :param path: path to csv file
    :return: dictionary of accession: ncbi_lineage
    """
    ncbi = NCBITaxa()
    d = {}
    with open(path) as fin:
        reader = csv.reader(fin, delimiter=sep)
        for row in reader:
            d[row[1]] = [str(x) for x in ncbi.get_lineage(row[0])]
    return d


def Counter(x):
    """
    Fast defaultdict Counter implementation
    Faster than collections.Counter
    """
    counter = defaultdict(int)
    for k in x:
        counter[k] += 1
    return counter


def percentage_sets(sets, percent=50, atmost=500):
    """
    function to return sets of values if present in a fraction of the sets

    :param sets: list of sets to reduce
    :param pecent: percentage of prevalence needed to include a value
    :param atmost: reduce final set to at most n values

    :return: set of values
    """
    n = len(sets)
    universe = []
    keep = set()
    for s in sets:
        universe.extend(list(s))
    C = Counter(universe)
    for element, count in C.items():
        if (count / n * 100) >= percent:
            keep.add(element)
    # reduce set to at most n elements
    # sorting and seting of the seed, will enforce reproducibility
    if len(keep) > atmost and atmost > 0:
        seed(125465)  # seed required to obtain the same SCMGs when reducing size, leading to reproducible results
        x = list(keep)
        x.sort()
        keep = set(sample(x, atmost))

    return keep


def union_sets(sets):
    """
    Given a list of sets will return the union between
    all sets

    :param sets: list of sets
    :type sets: lst

    :return: set with the union of all sets
    :rtype: set
    """
    if type(sets) is not list:
        raise TypeError("sets should be a list of sets")
    s = sets[0]
    for ss in sets:
        s = s.intersection(ss)
    return s


def load_SCMGs(path, sep=","):
    """
    Function to load and return all SCMGs found in a gzipped
    or not gipped file as a dictionary of sets
    with first column as key


    :param path: path to a two column csv file
    :type path: str
    :param sep: delimiter passed to `csv.reader`
    :type sep: str, optional


    :return: dictionary
    :rtype: dict
    """
    if not os.path.exists(path):
        raise FileNotFoundError("Could not find SCMG file")

    scmg = {}
    if path.endswith(".gz"):
        scmg_csv = gzip.open(path, mode="rt")
    else:
        scmg_csv = open(path)

    reader = csv.reader(scmg_csv, delimiter=sep)
    for row in reader:
        if len(row) != 2:
            raise ValueError("Expected exactly two columns")
        node, profile = row
        if node in scmg.keys():
            scmg[node].add(profile)
        else:
            scmg[node] = set([profile])

    scmg_csv.close()

    return scmg


def which(program):
    """
    test if w programm is avaliable

    :param programm: name of an executable

    :rtype: bool
    """
    # taken from
    # https://stackoverflow.com/questions/377017/\
    # test-if-executable-exists-in-python
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


def read_hmmer(path, cutoffs_p=None):
    """
    Read in a hmmer output file and return dictionary of columns
    """
    if not os.path.exists(path):
        raise FileNotFoundError("Could not find hmmer file")

    # load cutoffs from csv if avail
    cutoffs = defaultdict(str)
    if cutoffs_p is not None:
        with open(cutoffs_p) as fin:
            for row in csv.DictReader(fin):
                cutoffs[row["query"]] = row["GA"]

    columns = {"target": 0, "query": 2, "bitscore": 5, "evalue": 4}
    data = []
    with open(path) as hin:
        for line in hin:
            if line.startswith("#"):
                continue
            row = line.strip().split()
            d = {k: row[v] for k, v in columns.items()}
            if cutoffs_p is not None:
                d["expected_GA"] = cutoffs[row[columns["query"]]]
            data.append(d)

    return data


def hmmer_cutoffs(path, outfile):
    """
    Read in a hmmer input file and return dictionary of columns
    """
    if not os.path.exists(path):
        raise FileNotFoundError("Could not find hmm file")

    name = False
    with open(path) as hin, open(outfile, "w") as fout:
        of = csv.DictWriter(fout, delimiter=",", fieldnames=["query", "GA"])
        of.writeheader()
        # parse hmm file
        for line in hin:
            if line.startswith("#"):
                continue
            l = line.split()
            if l[0] in ["NAME", "GA"]:
                if l[0] == "NAME":
                    name = l[1]
                elif name is not False:
                    of.writerow({"query": name, "GA": l[1]})
    return outfile


def horizontal_concatenate(output, files, profiles):
    """
    Horizontally concatenate fasta files

    """

    def read_fasta(path):
        seqs = {}
        for record in Fasta(path):
            seqs[record.name] = record.seq
        return seqs

    # make sure profiles are sorted
    profiles.sort()
    seqs = {}
    # first we detect the set of all names
    allnames = set()
    # determine the length of each alignment for padding
    lengths = defaultdict(int)
    for profile, f in zip(profiles, files):
        seqs[profile] = read_fasta(f)
        for name, seq in seqs[profile].items():
            allnames.add(name)
            # remember the lengths of this profile
            lengths[profile] = len(seq)

    # add missing seqences as strings of spaces
    for profile in seqs.keys():
        for name in allnames:
            if name not in seqs[profile].keys():
                seqs[profile][name] = "".join(lengths[profile] * ["-"])

    mergedSeqs = {}
    for name in allnames:
        s = ""
        for p in profiles:
            s += seqs[p][name]
        mergedSeqs[name] = s

    # write output:
    with open(output, "w") as f:
        for name, seq in mergedSeqs.items():
            seq = seq.replace(".", "-")
            f.write(">{}\n{}\n".format(name, seq))

    return output
