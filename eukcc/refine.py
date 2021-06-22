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


def connected_bins(name, link_table, min_links):
    connected = set()
    for name1, d in link_table.items():
        for name2, links in d.items():
            if name1 == name2:
                continue
            if name1 == name or name2 == name:
                if links >= min_links:
                    n = name2 if name1 == name else name1
                    connected.add(n)
    return connected


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
        logging.debug(
            "Using existing merged contigs, delete output folder if thats not what you want."
        )

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
    # s_cmb = n_combi(smallbins, state["n_combine"])
    # logging.info("Created {} possible combinations of small bins".format(len(s_cmb)))

    n_large_bins = len(bins) - len(smallbins)
    logging.info("Found {} large bins to merge with".format(n_large_bins))

    refined = []

    def get_name(i):
        return bins[i].name.split("_", 1)[1]

    def combine_bins(idx, children, smallbins, stop=0):
        all_comb = []

        possible_kids = []
        if state["links"] is not None:
            possible_kids = []
            connected_kids = connected_bins(
                get_name(idx), link_table, state["min_links"]
            )
            for k in children:
                connected_kids = connected_kids.union(
                    connected_bins(get_name(k), link_table, state["min_links"])
                )
            # convert names back to indicies
            for i, b in enumerate(bins):
                if get_name(i) in connected_kids:
                    possible_kids.append(i)
        else:
            possible_kids = smallbins

        possible_kids = set(possible_kids) - set(children)

        # add a single kid and evaluate again
        if stop < 1:
            s = list(set(children) | set([idx]))
            s.sort()
            return [tuple(s)]
        else:
            for kid in possible_kids:
                k = children.copy()
                k.append(kid)
                all_comb.extend(combine_bins(idx, k, smallbins, stop=stop - 1))
            return all_comb

    # mergers = defaultdict(lambda: {'parent': None, "children": [], "gain": []})
    parent_bins = [i for i, b in enumerate(bins) if i not in smallbins]
    for parent_bin in parent_bins:
        parent_name = bins[parent_bin].name.split("_", 1)[1]

        # we find all combinations that are valid to merge with
        combies = []
        for i in range(state["n_combine"]):
            ks = combine_bins(parent_bin, [], smallbins, stop=i)
            combies.extend(ks)
        combies = list(set(combies))
        # sort combinations by length
        combies = sorted(combies, key=len, reverse=False)

        logging.info(
            "For bin {} we found {} merging combinations".format(
                parent_name, len(combies)
            )
        )
        for kid_ids in combies:
            # turn kid indexec into sorted list, so they are reproducible turn kid indexec into sorted list, so they are reproducible
            kid_ids = list(kid_ids)
            kid_ids.sort()
            kid_ids = [i for i in kid_ids if i != parent_bin]
            children = [bins[i] for i in kid_ids]
            if len(children) == 0:
                continue
            logging.debug("Testing bin combination for bin {}".format(parent_name))
            logging.debug(
                "Adding in bins {}".format(",".join([get_name(i) for i in kid_ids]))
            )
            # make all merges and see if they are great or not
            merged = merge_bins(
                bins[parent_bin], children, os.path.join(state["workdir"], "mergers")
            )

            # identify parent bin, if [A, B, C] then parent can be [A,B] or [A,C] or if none of them exists
            # it should be [A]
            compare_state = bins[parent_bin].state
            if len(children) > 1:
                for c in combinations(kid_ids, len(kid_ids) - 1):
                    for potential_parent in refined:
                        if set(c) == set(potential_parent["children_idx"]):
                            compare_state = potential_parent
                            break

            gain_cp = round(
                merged["quality"]["completeness"]
                - compare_state["quality"]["completeness"],
                2,
            )
            gain_ct = round(
                merged["quality"]["contamination"]
                - compare_state["quality"]["contamination"],
                2,
            )
            if gain_cp >= state["improve_percent"] and gain_cp > (
                state["improve_ratio"] * gain_ct
            ):
                logging.warning(
                    "Successfull merge: +{}% Compl. +{}% Cont.".format(gain_cp, gain_ct)
                )
                merged["children_idx"] = kid_ids
                merged["parent_idx"] = parent_bin
                refined.append(merged)

    # logging.info("Iterated all possible combinations")
    # refined = remove_double_kids(refined, bins)

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
