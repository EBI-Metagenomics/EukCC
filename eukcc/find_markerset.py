import os
import argparse
import logging
from collections import defaultdict, Counter
from eukcc.eukcc import eukcc, eukcc_state
import eukcc.version as version
from eukcc.file import file
import math
from multiprocessing import Pool
from functools import partial
from operator import itemgetter
from random import sample


def fake(fasta, state):
    return fasta


def search_genome(fasta, state):
    """
    search a genome for all markers in a clade
    """
    logging.info("Running for {}".format(fasta))

    workdir = os.path.join(state["workdir"], os.path.basename(fasta))
    my_state = eukcc_state(workdir, state)
    my_state["fasta"] = fasta
    E = eukcc(my_state)
    E.predict_protein()
    markers = E.search_markers(use_all=True)
    logging.info("Done")
    return markers


def score_profile(lst, ignore=None, lenient=False):
    """
    score a profile based on how many SCMG are contained
    if many score is high, if few, score is low
    we can also pass a number of accesions to ignore for the score

    The second score traces how many singletons in total would be covered
    This is to rank profiles even if the number of novel singleton coverage
    is the same, the hidden coverage (Of genomes that are already
    well covered) might be different
    """
    if ignore is None:
        ignore = set()
    # in lenient mode we allow also doubletons to be counted as a win
    if lenient:
        max_cov = 2
    else:
        max_cov = 1
    score = 0
    hidden_score = 0
    for key, value in Counter(lst).items():

        increase = 0
        if value <= max_cov:
            increase += 1
            if value == 1:
                hidden_score += 1
        elif value > 2:
            # doubletons get no negative penalty
            increase += -value

        # if we have enough markers for a given species,
        # any marker gets a score of 0
        if key in ignore and increase > 0:
            increase = 0
        score = score + increase
    return (score, hidden_score)


def greedy_marker_selection(
    all_profiles, ignore_profiles=None, ignore_genomes=None, max_iter=50, target_cov=20
):
    if ignore_genomes is None:
        ignore_genomes = []
    if ignore_profiles is None:
        ignore_profiles = []
        selected_profile = []
    else:
        selected_profile = ignore_profiles.copy()

    # list of all genomes
    universe = set()
    # list of all profiles
    profile_universe = set()
    for profile, hits in all_profiles.items():
        universe = universe | set(hits)
        if profile not in ignore_profiles:
            profile_universe.add(profile)
    # remove ignored genomes
    universe = universe - set(ignore_genomes)
    # transform to list, so its stable in itration
    profile_universe = list(profile_universe)

    # selected profiles in this round
    covered_genomes = set()
    covered_genomes = covered_genomes | set(ignore_genomes)
    score = defaultdict(int)

    for i in range(0, max_iter):
        logging.debug(
            "Found enough profiles for {} out of {}".format(
                len(covered_genomes), len(universe)
            )
        )
        if len(covered_genomes) >= len(universe):
            logging.info("Found enough profiles to cover everything")
            break

        # score all profiles and sort them
        scores = [
            score_profile(all_profiles[profile], covered_genomes)
            for profile in profile_universe
        ]

        scores_sort = sorted(scores, key=itemgetter(0, 1), reverse=True)
        max_idx = [idx for idx, scores in enumerate(scores) if scores == scores_sort[0]]
        idx = sample(max_idx, 1)[0]

        if scores_sort[0][0] < 0.05 * (len(universe) - len(covered_genomes)):
            logging.info(
                "No profile covers more than 5% of the remaining universe, so we stop here"
            )
            break

        if scores_sort[0][1] < 5:
            logging.warning("No profile covers enough (5) genomes anymore")
            i = 0
            for s in scores_sort:
                print(s)
                i += 1
                if i > 5:
                    break
            break

        best_profile = profile_universe[idx]
        selected_profile.append(best_profile)
        del profile_universe[idx]

        # make sure we can ignore covered genomes in scoring
        for acc in set(all_profiles[best_profile]):
            score[acc] += 1
        for key, value in score.items():
            if value >= target_cov:
                covered_genomes.add(key)

    # compute coverage stats from last score
    stats = {}
    for acc, cov in score.items():
        # identify multitons
        n = 0
        for profile in selected_profile:
            if all_profiles[profile].count(acc) > 1:
                n += 1
        logging.debug(
            "Covered {} with {} profiles, {} of which are multitons".format(acc, cov, n)
        )
        stats[acc] = {"covered": cov, "multitons": n}
    return (selected_profile, stats)


def define_tree_set(data, n_target=30):
    all_profiles = defaultdict(list)
    universe = set()
    for i, d in enumerate(data):
        universe.add(i)
        for row in d:
            all_profiles[row["query"]].append(i)

    # first round
    choosen_profiles, stats = greedy_marker_selection(
        all_profiles,
        ignore_profiles=None,
        ignore_genomes=None,
        max_iter=n_target * 3,
        target_cov=round(n_target * 0.7),
    )

    ignore_genomes = set()
    for acc, stat in stats.items():
        if stat["covered"] >= round(n_target * 0.7):
            ignore_genomes.add(acc)
    # second round
    choosen_profiles, stats = greedy_marker_selection(
        all_profiles,
        ignore_profiles=choosen_profiles,
        ignore_genomes=ignore_genomes,
        max_iter=n_target,
        target_cov=round(n_target * 0.5),
    )

    ignore_genomes = set()
    for acc, stat in stats.items():
        if stat["covered"] >= n_target:
            ignore_genomes.add(acc)
    # final_round
    choosen_profiles, stats = greedy_marker_selection(
        all_profiles,
        ignore_profiles=choosen_profiles,
        ignore_genomes=ignore_genomes,
        max_iter=n_target,
        target_cov=round(n_target / 2),
    )
    logging.info(
        "Defined {} profiles for tree construction".format(len(choosen_profiles))
    )

    return choosen_profiles


def find_intersection(data, missing=50):
    singletons = defaultdict(int)
    multitons = defaultdict(int)

    for d in data:
        found = Counter([x["query"] for x in d])
        for profile, number in found.items():
            multitons[profile] += 1
            if number == 1:
                singletons[profile] += 1
    # construct a result
    results = defaultdict(set)
    for i in range(0, missing):
        logging.debug("constructing set for {} missing".format(i))
        n_genomes = len(data) - i
        for profile, count in singletons.items():
            if count == n_genomes:
                results[n_genomes].add(profile)
    return results


def main():
    # set arguments
    # arguments are passed to classes
    parser = argparse.ArgumentParser(
        description="Evaluate completeness and contamination of a MAG."
    )
    parser.add_argument(
        "genomes", type=str, help="Find marker for these genomes", nargs="+"
    )
    parser.add_argument(
        "--out",
        "-o",
        type=str,
        required=False,
        help="Path to output folder (Default: .)",
        default=".",
    )
    parser.add_argument("--db", type=str, default=None, help="Path to EukCC DB")
    parser.add_argument(
        "--threads", type=int, help="Number of threads to use (Default: 1)", default=1
    )
    parser.add_argument(
        "--tree",
        type=int,
        help="Number of profiles to use at target for tree profiles (default: 30)",
        default=30,
    )
    parser.add_argument(
        "--clade",
        default="base",
        type=str,
        help="Define clade as base, fungi, protozoa or plants (Defaut: base)",
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
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="EukCC version {}".format(version.__version__),
    )
    args = parser.parse_args()
    state = eukcc_state(
        workdir=os.path.join(args.out, "refine_workdir"), options=vars(args)
    )
    file.isdir(state["workdir"])

    # define logging
    logLevel = logging.INFO
    if state["quiet"]:
        logLevel = logging.WARNING
    elif state["debug"]:
        logLevel = logging.DEBUG
    logging.basicConfig(
        format="%(asctime)s %(message)s",
        datefmt="%d-%m-%Y %H:%M:%S: ",
        level=logLevel,
    )
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

    logging.info("EukCC version {}".format(version.__version__))

    logging.info(
        "Looking for shared markers across {} genomes".format(len(state["genomes"]))
    )
    n_per_worker = 4  # using more threads for hmmer makes no sense, so we parallize accroos genomes
    if state["threads"] > (2 * n_per_worker):
        # multithreading pool
        n_processes = math.floor(state["threads"] / n_per_worker)
        logging.info(
            "Launching {} threads with {} threads each".format(
                n_processes, n_per_worker
            )
        )
        pool = Pool(processes=n_processes)
        # change threads not
        opt = {k: v for k, v in state.opt.items()}
        opt["threads"] = n_per_worker
        search_genome_p = partial(search_genome, state=opt)
        data = pool.map(search_genome_p, state["genomes"])
        pool.close()
        pool.join()
    else:
        data = []
        for genome in state["genomes"]:
            data.append(search_genome(genome, state))

    tree_profiles = define_tree_set(data, n_target=args.tree)
    result = find_intersection(data, missing=3)
    outfile = os.path.join(state["out"], "profiles.txt")
    with open(outfile, "w") as fout:
        for key, profiles in result.items():
            for profile in profiles:
                fout.write("{}\t{}\n".format(key, profile))
        for profile in tree_profiles:
            fout.write("{}\t{}\n".format("tree", profile))
    logging.info("wrote profiles to {}".format(outfile))
