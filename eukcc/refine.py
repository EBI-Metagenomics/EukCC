import os
import logging
import glob
import csv
from collections import defaultdict
from eukcc.eukcc import eukcc, eukcc_state
from eukcc.fasta import merge_fasta
import eukcc.version as version
from eukcc.fasta import Fasta
from eukcc.bin import bin, merge_bins
from eukcc.file import file
from itertools import combinations


def split_contig_faa(path, workdir, delim="_binsep_"):
    faa_dir = os.path.join(workdir, "faa")
    if os.path.exists(faa_dir):
        logging.info(
            "Faa folder exists, remove if you want to rerun this analsys, will reuse it for now"
        )
        return [os.path.abspath(os.path.join(faa_dir, x)) for x in os.listdir(faa_dir)]
    elif file.isdir(faa_dir):
        fls = {}
        for seq in Fasta(path):
            b = seq.name.split(delim, 1)[0]
            contig = seq.name.split(delim, 1)[1]
            if b not in fls.keys():
                fls[b] = open(os.path.join(faa_dir, b), "w")
            fls[b].write(">{name}\n{seq}\n".format(name=contig, seq=seq.seq))
        # close all files
        for key, fl in fls.items():
            fl.close()
        return [os.path.abspath(os.path.join(faa_dir, x)) for x in os.listdir(faa_dir)]
    else:
        logging.warning("Could not create fasta split dir: {}".format(faa_dir))


def eukcc_folder(args):
    state = eukcc_state(
        workdir=os.path.join(args.out, "refine_workdir"), options=vars(args)
    )
    file.isdir(state["workdir"])

    if state["db"] is None:
        if os.environ.get("EUKCC2_DB") is not None:
            state["db"] = os.environ.get("EUKCC2_DB")
            logging.debug(
                "Defined db via env variable EUKCC2_DB as '{}'".format(state["db"])
            )
        else:
            logging.error("No database was provided via --db or EUKCC2_DB env variable")
            exit(202)

    logging.info("EukCC version {}".format(version.__version__))

    # discover all bins
    logging.debug("Loading all bins")
    state["input_bins"] = glob.glob(
        os.path.join(state["binfolder"], "*{}".format(state["suffix"]))
    )
    logging.info("Found {} bins".format(len(state["input_bins"])))
    if len(state["input_bins"]) == 0:
        logging.error("Stopping as no bins in folder")
        exit(222)
    refine(state)


def check_links(names, link_table, min_links=100):
    for name_1 in names:
        for name_2 in names:
            if name_1 == name_2:
                continue
            else:
                if link_table[name_1][name_2] < min_links:
                    return False
    return True


def read_link_table(path):
    bin_table = defaultdict(lambda: defaultdict(int))
    with open(path) as fin:
        for row in csv.DictReader(fin):
            bin_table[row["bin_1"]][row["bin_2"]] = int(row["links"])
            bin_table[row["bin_2"]][row["bin_1"]] = int(row["links"])
    return bin_table


def write_table(state, path, sep="\t", header=False):
    if header:
        with open(path, "w") as fout:
            fout.write(sep.join(["bin", "completeness", "contamination"]))
            fout.write("\n")
    else:
        if state["marker_set"]["quality"] != "good":
            return
        with open(path, "a") as fout:
            bn = os.path.basename(state["fasta"])
            # trim metaeuk prefix
            if bn.startswith("metaeuk_"):
                bn = bn.split("_", 1)[1]
            compl = str(state["quality"]["completeness"])
            cont = str(state["quality"]["contamination"])
            s = sep.join([bn, compl, cont])
            fout.write(s)
            fout.write("\n")


def refine(state):
    """
    main refinement pipeline
    It takes a state argument and does the rest. It is basically a per bin EukCC run
    with added mergin features
    """
    state["contigs"] = os.path.join(state["workdir"], "contigs.fasta")
    if not os.path.exists(state["contigs"]):
        logging.debug("Merging bins into temp contig file")
        state["contigs"] = merge_fasta(
            state["input_bins"], state["contigs"], seperator="_binsplit_"
        )
    else:
        logging.debug("Using existing merged contigs")

    # initialize output tables
    result_table = os.path.join(state["out"], "eukcc.csv")
    write_table(None, result_table, header=True)

    merged_table = os.path.join(state["out"], "merged_bins.csv")
    note_merges(None, None, merged_table, header=True)

    # predict proteins using metaeuk
    state["fasta"] = state["contigs"]
    logging.debug("Initialize EukCC")
    E = eukcc(state)
    E.predict_protein()
    # split proteins into Bin files
    # this also gets rid of bins with zero proteins
    state["faas"] = split_contig_faa(state["faa"], state["workdir"], delim="_binsplit_")

    if state["links"] is not None:
        link_table = read_link_table(state["links"])

    # for each bin create a EukCC run and try to place and estimate its completeness
    bins = []
    for i, path in enumerate(state["faas"]):
        logging.debug(
            "Running EukCC on: {} ({}/{})".format(
                os.path.basename(path), i, len(state["faas"])
            )
        )
        wd = os.path.join(
            state["workdir"], "refine", "bin_{}".format(os.path.basename(path))
        )
        bins.append(bin(state, wd, path, protein=True))

    smallbins = []
    for i, b in enumerate(bins):
        if b.state["quality"] is None:
            smallbins.append(i)
            continue
        elif b.state["quality"]["completeness"] < 50:
            smallbins.append(i)
            write_table(b.state, result_table)
            continue
        else:
            write_table(b.state, result_table)

    # get all combinations of small bins
    s_cmb = n_combi(smallbins, state["n_combine"])
    logging.info("Created {} possible combinations of small bins".format(len(s_cmb)))

    n_large_bins = len(bins) - len(smallbins)
    logging.info("Found {} large bins to merge with".format(n_large_bins))

    max_iters = len(bins)
    already_at = 0
    refined = []
    for i, b in enumerate(bins):
        already_at += 1
        logging.debug(
            "Refining bin {}/{} {}%".format(
                already_at, max_iters, round(100 * already_at / max_iters, 1)
            )
        )
        if i in smallbins:
            continue
        if (
            bins[i].state["quality"]["completeness"] < 50
            or bins[i].state["quality"]["contamination"] >= 15
            or bins[i].state["marker_set"]["quality"] != "good"
        ):
            # skip non medium quality bins
            continue
        # big bins can be combined with a set of small bins
        # each big bin we will merge the faa with the small bin faas
        # then we use the set of the big bin to estimate the merger
        # we retain the merger if it fits the layed out parameters
        for s in s_cmb:
            # construct children container
            children, names, skip = build_children(bins, s, i)
            if skip:
                continue
            # merge bins
            logging.debug("Testing combination {}".format(names))
            # check for linking reads
            if state["links"] is not None:
                if check_links(names, link_table, state["min_links"]) is False:
                    logging.debug("Skipping merger with to few supporting reads")
                    continue

            merged = merge_bins(
                bins[i], children, os.path.join(state["workdir"], "mergers")
            )
            gain_cp = round(
                merged["quality"]["completeness"]
                - bins[i].state["quality"]["completeness"],
                2,
            )
            gain_ct = round(
                merged["quality"]["contamination"]
                - bins[i].state["quality"]["contamination"],
                2,
            )
            if gain_cp >= state["improve_percent"] and gain_cp > (
                state["improve_ratio"] * gain_ct
            ):
                logging.debug(
                    "Successfull merge: +{}% Compl. +{}% Cont.".format(gain_cp, gain_ct)
                )
                merged["children_idx"] = s
                merged["parent_idx"] = i
                refined.append(merged)

    logging.info("Iterated all possible combinations")
    refined = remove_double_kids(refined, bins)

    # Here we just have to export the bins at this point. Thats easy and quick.
    # For each new bin we rerun eukcc with a enw metaeuk run, so we get a final verdict
    ref_dir = os.path.join(state["out"], "merged_bins")
    if os.path.exists(ref_dir):
        logging.warning(
            "Folder merged_bins alerady exist. Will not overwrite. Please remove the folder."
        )
        merged_fnas = [os.path.join(ref_dir, x) for x in os.listdir(ref_dir)]
    else:
        file.isdir(ref_dir)
        merged_fnas = []
        for merged_i, merged in enumerate(refined):
            new_name = "{}{}{}".format(state["prefix"], merged_i, state["suffix"])
            merged_fna = os.path.join(ref_dir, new_name)
            # combined index list
            idxs = [merged["parent_idx"]]
            for i in merged["children_idx"]:
                idxs.append(i)
            names = [os.path.basename(bins[i].state["faa"]) for i in idxs]
            names = [n.split("_", 1)[1] for n in names]
            # write that to file
            note_merges(new_name, names, merged_table)

            fnas = [os.path.join(state["binfolder"], name) for name in names]
            merged_fna = merge_fasta(fnas, merged_fna, seperator=None)
            merged_fnas.append(merged_fna)
            logging.info("Created combined fasta {}".format(merged_fna))

    evaluate_multiples(merged_fnas, result_table, state)
    logging.info("Created {} merged bins".format(len(refined)))


def note_merges(new_name, names, path, header=False, sep="\t", child_sep="child:"):
    if header:
        with open(path, "w") as fout:
            fout.write(sep.join(["merged", "bins"]))
            fout.write("\n")
    else:
        with open(path, "a") as fout:
            children = child_sep.join(names)
            s = sep.join([new_name, children])
            fout.write(s)
            fout.write("\n")


def build_children(bins, subset, i):
    skip = False
    children = []
    for indx in subset:
        if bins[indx].state["clade"] is not None:
            bin_clade = bins[indx].state["clade"]
            parent_clade = bins[i].state["clade"]
            if (
                parent_clade != "base" and bin_clade != "base"
            ) and bin_clade != parent_clade:
                # parent and child have diverging clades, not need to merge any of these
                logging.debug(
                    "Parent and child clade diverge: {}/{}".format(
                        parent_clade, bin_clade
                    )
                )
                skip = True
        children.append(bins[indx])

    # construct names
    names = []
    names.append(bins[i].name)
    for c in children:
        names.append(c.name)
    names = [n.split("_", 1)[1] for n in names]

    return children, names, skip


def evaluate_multiples(fnas, result_table, state):
    # launch eukcc on all bins
    for path in fnas:
        wd = os.path.join(
            state["workdir"], "refine", "bin_{}".format(os.path.basename(path))
        )
        E = bin(state, wd, fasta=path, protein=False)
        write_table(E.state, result_table)


def remove_double_kids(refined, bins):
    ####################################################################
    # Initially make sure each child is used only for a single parent
    child_parents = defaultdict(list)
    remove_child = set()
    for merged in refined:
        for child in merged["children_idx"]:
            child_parents[child].append(merged["parent_idx"])
    for child, parents in child_parents.items():
        if len(set(parents)) > 1:
            logging.info(
                "Child bin {} as assigned multiple parents. Thus they will be ignored to avoid hybrids".format(
                    bins[child].state["name"]
                )
            )
            remove_child.add(child)
    del_idx = [
        indx
        for indx, merged in enumerate(refined)
        if len(remove_child & set(merged["children_idx"])) > 0
    ]
    for indx in del_idx:
        del refined[indx]
        for indx in del_idx:
            del refined[indx]
    ####################################################################
    return refined


def n_combi(smallbins, n):
    s_cmb = []
    for i in range(1, n + 1):
        s_cmb.extend(list(combinations(smallbins, i)))
    s_cmb = list(set(s_cmb))
    s_cmb.sort()
    return s_cmb
