import os
import logging
from collections import defaultdict


def check_alignment(path, target, target_upper=10):
    upper = 0
    for seq in Fasta(target):
        target_id = seq.name
        break

    for seq in Fasta(path):
        if seq.name != target_id:
            continue
        upper = sum(1 for c in seq.seq if c.isupper())
    return upper >= target_upper


def determine_type(path):
    """
    look at the first N entries to figuer out if its a DNA entry
    """
    n = 5
    seqtype = "DNA"
    valid_DNA_chars = set(list("AGCTURYKMSWBDHVN"))
    fa = Fasta(path)
    for record in fa:
        if len(set(list(record.seq.upper())) - valid_DNA_chars) > 0:
            # this seems to not be a DNA sequence
            seqtype = "AA"
            break

        if n < 0:
            break
        else:
            n = n - 1

    return seqtype


def merge_fasta(files, outfile, seperator=None):
    """
    merges faster and adds a seperator if wanted
    """
    with open(outfile, "w") as fout:
        for fasta in files:
            binname = os.path.basename(fasta)
            for seq in Fasta(fasta):
                if seperator is not None:
                    name = "{}{}{}".format(binname, seperator, seq.name)
                else:
                    name = seq.name
                fout.write(">{name}\n{seq}\n".format(name=name, seq=seq.seq))
    return outfile


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


def N50(path):
    """
    Will compute N50 for a given Fasta file
    """
    sizes = []
    for rec in Fasta(path):
        sizes.append(len(rec.seq))
    sizes.sort(reverse=True)
    complete_length = sum(sizes)
    rs = 0
    for s in sizes:
        rs = rs + s
        if rs >= complete_length / 2:
            break
    return s


def Fasta(path):
    """
    Iterator for fasta files
    """
    entry = False
    if validate_fasta(path) is False:
        raise ValueError("The provided fasta file is malformed: {}".format(path))

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


def reduce_fasta(path, outfile, ids):
    if type(ids) is not list:
        raise TypeError("Ids must be of list type")

    with open(outfile, "w") as fout:
        for seq in Fasta(path):
            if seq.name in ids or seq.long_name in ids:
                fout.write(">{}\n{}\n".format(seq.long_name, seq.seq))
    return outfile


def validate_fasta(path, seqtype="AA"):
    """
    function to make sure fasta files are corretly formated
    not empty and not a product of a failed GeneMark-ES attempt

    :params path: fasta file
    :type path: str

    :rtype: bool
    """
    if not os.path.exists(path):
        raise FileNotFoundError("Could not find fasta file")

    # empty files are not valid
    if os.stat(path).st_size == 0:
        return False

    # check if we can open it
    try:
        with open(path) as fin:
            for line in fin:
                break
    except UnicodeDecodeError:
        return False

    # check if the file contains valid sequences
    with open(path) as fin:
        tp = "start"
        for line in fin:
            if line.strip() == "":
                continue
            if tp == "start" and not line.strip().startswith(">"):
                # the first non empty line needs to be a name
                logging.error("Fasta file contains non name line as first non empty line")
                return False

            if line.strip().startswith(">") and tp == "name":
                # we can not have two names in a row
                logging.error("Fasta contains two names in row")
                return False
            elif line.strip().startswith(">"):
                # if we had sequence before we can have a name now
                tp = "name"

            if not line.strip().startswith(">"):
                tp = "seq"

    return True


def clean_metaeuk_fasta(fasta, outfasta):
    """
    given a fasta file, this will cleanup names, so we dont have any spaces or
    characters that confuse epa-ng such as  |[]
    """

    nms = defaultdict(int)
    with open(outfasta, "w") as fout:
        # in case of an empty fasta, we also return an emopty fasta
        if os.stat(fasta).st_size == 0:
            return outfasta
        for seq in Fasta(fasta):
            # clean name using contig as the name
            seq.name = seq.name.split("|")[1]
            # we use name_ORF as name format, this could be compatible with CAT
            name = "{name}_{i}".format(name=seq.name, i=nms[seq.name])
            # increase name counter
            nms[seq.name] += 1
            fout.write(">metaeuk_{name}\n{seq}\n".format(name=name, seq=seq.seq))
    return outfasta
