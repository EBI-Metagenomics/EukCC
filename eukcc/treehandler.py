from ete3 import Tree
from collections import defaultdict, Counter
import logging
from eukcc.base import load_SCMGs, percentage_sets, load_tax_info


class markerset:
    def __init__(self, leafes, covered, profiles, prevalence=100):
        self.leafes = leafes
        self.profiles = profiles
        self.covered = covered
        self.nspecies = len(leafes)
        self.nprofiles = len(profiles)
        self.prevalence = prevalence
        # self.all_leafes = None
        self.all_places = None
        self.density = len(self.covered) / self.nspecies

    def to_dict(self):
        # estimate quality for easy usage
        self.quality = "good"
        if len(self.covered) < (0.5 * len(self.all_places)):
            self.quality = "bad"
        return {
            "profiles": self.profiles,
            "prevalence": self.prevalence,
            "covered": len(self.covered),
            "placements": len(self.all_places),
            "nspecies": self.nspecies,
            "quality": self.quality,
        }

    @property
    def score(self):
        # initial model, with very little training data for comaprison
        # pred_error = 26.465607 + -0.252184 * self.prevalence + 0.004341 * len(self.leafes) + 0.008279 * len(self.profiles) + -12.807798 * (len(self.covered)/len(self.all_places))
        # pred_error = 36.131638 + -0.263107 * self.prevalence + -0.011279 * len(self.leafes) + -0.001418 * len(self.profiles) + -13.084319 * (len(self.covered)/len(self.all_places))
        pred_error = (
            133.524566
            + -0.458508 * self.prevalence
            + -0.038799 * len(self.leafes)
            + -0.065516 * len(self.profiles)
            + -35.417479 * (len(self.covered) / len(self.all_places))
        )
        # the biggest score has the lowest absolute error
        score = abs(pred_error)

        return score

    @property
    def best_guess(self):
        """
        Handcrafted priorty list to get a  good marker gene set
        """
        if len(self.all_places) < 10:
            cover_score = 10 * len(self.covered)
        else:
            cover_score = 100 * (len(self.covered) / len(self.all_places))

        if self.nprofiles > 200:
            profile_score = 100
        elif self.nprofiles > 100:
            profile_score = 80
        elif self.nprofiles > 30:
            profile_score = 60
        else:
            profile_score = self.nprofiles

        prevalence_score = (self.prevalence - 100) * 10

        ns_t = 7
        if self.nspecies > ns_t:
            species_score = 100
        else:
            species_score = 100 * (self.nspecies / ns_t)

        return cover_score + profile_score + prevalence_score + species_score

    def __str__(self):
        return "set with score {} and {} profiles {}/{}/{}".format(
            self.score,
            len(self.profiles),
            self.prevalence,
            len(self.covered),
            len(self.leafes),
        )


def hard_set_computation(set_path, genomes, prevalence=98, atmost=500, set_size=20):
    """
    Function to compute set based on a list of genomes passed to it
    """
    scmg = load_SCMGs(set_path)
    found = False
    set_prevalence = 100
    biggest = 0
    while found is False and set_prevalence >= prevalence:
        logging.debug(
            "Searching for Marker set at {} prevalence across {} genomes".format(
                set_prevalence, len(genomes)
            )
        )
        sets = []
        for genome in genomes:
            try:
                sets.append(scmg[genome])
            except KeyError:
                logging.warning(
                    "Database missing markes for '{}'. This should not be the case. Make sure the database is not corrupted".format(
                        genome
                    )
                )
        s = percentage_sets(sets, set_prevalence, atmost)
        if len(s) > biggest:
            biggest = len(s)
        if len(s) >= set_size:
            found = True
            break
        set_prevalence = set_prevalence - 0.5

    logging.debug("Largest set we found had {} SCMGs".format(biggest))
    if found:
        logging.debug(
            "Found set of size {} with prevalence {}".format(len(s), set_prevalence)
        )
        return s
    else:
        return None


def tax_LCA(tree, taxinfo, placements=None, majority_vote=0.6):
    """
    given a guppy tree we will try to figure out
    the placements
    """
    logging.debug("Finding LCA")
    t = Tree(tree)
    info = load_tax_info(taxinfo)

    if placements is None:
        known_leafes = set(info.keys())
    else:
        node = t.get_tree_root()
        known_leafes = set(node.get_leaf_names()) - set(placements)

    logging.debug(
        "Determining the LCA in a tree with {} leafes".format(len(known_leafes))
    )

    leafes = set()
    all_leafes = set(node.get_leaf_names())
    pl_lngs = []
    for place in all_leafes - known_leafes:
        # jump to node
        node = t.search_nodes(name=place)[0]
        while True:
            # reduce to possible keys, removing placements
            leafes = set(node.get_leaf_names()) & known_leafes
            if len(leafes) > 0:
                break
            else:
                node = node.up
        lngs = [info[acc].copy() for acc in leafes]

        def maj_lng(L, f=0.6, add_protists=True):
            non_protists = ["33154", "33090"]
            lngs = L.copy()
            # get the tax ids at each position
            positions = defaultdict(list)
            index = set()
            for lng in lngs:
                for i, tax in enumerate(lng):
                    if i == 3 and tax not in non_protists and add_protists:
                        # insert protists shared taxid
                        positions[i].append("protist_common")
                    positions[i].append(tax)
                    index.add(i)
            idx = 0
            frac = 0
            lng = []
            while idx < max(index):
                C = Counter(positions[idx])
                idx += 1

                max_tax = C.most_common()[0][0]
                n = C[max_tax]
                frac = n / len(lngs)
                if frac >= majority_vote:
                    lng.append(max_tax)
                else:
                    logging.debug("majority vote to low")
                    break
            return lng

        pl_lngs.append(maj_lng(lngs))

    return maj_lng(pl_lngs)


class tree_sets:
    """
    Handle all tree operations

    :param path: path to newick tree file
    :type path: str

    """

    def __init__(
        self,
        tree_v,
        placement,
        setp,
        set_species=5,
        set_size=50,
        set_prevalence=98,
        set_atmost=500,
        dynamic_root=False,
        set_selection="lm",
        use_ncbi=False,
        training=False,
        taxinfo=None,
    ):
        self.t = Tree(tree_v)

        # find LCA of all placements
        # make places into convenient list
        pl = [x["n"] for x in placement["placements"]]
        places = []
        for p in pl:
            if type(p) is list:
                places.extend(p)
            elif type(p) is str:
                places.append(p)
        # root the tree to get best clade patterns
        if dynamic_root:
            logging.debug("Will use most distant entry to LCA as outgroup")
            self.lca = self.LCA(places)
            new_root = self.lca.get_farthest_node()
            self.t.set_outgroup(new_root[0])

        # load in all marker genes
        scmg = load_SCMGs(setp)
        self.known_leafes = set(load_tax_info(taxinfo).keys())

        logging.debug(
            "Starting to look for scmg set, selection based on {}".format(set_selection)
        )
        if use_ncbi:
            logging.debug("Will use NCBI tree instead of eukcc tree")
            self.marker_set = self._find_best_ncbi_set(
                places,
                scmg,
                taxinfo=taxinfo,
                min_set_size=set_size,
                set_atmost=set_atmost,
                set_species=set_species,
                min_prevalence=set_prevalence,
            )
        else:
            # expose final prevalence
            self.marker_set = self._find_best_set(
                places,
                scmg,
                training=training,
                min_set_size=set_size,
                set_atmost=set_atmost,
                set_species=set_species,
                min_prevalence=set_prevalence,
                sort_using=set_selection,
            )

        if training is False and self.marker_set is not None:
            logging.debug(
                "Defined SCMG set with {} marker genes with a single copy prevalence of {} percent covering {} related genomes supported by {}/{} placements".format(
                    len(self.marker_set.profiles),
                    self.marker_set.prevalence,
                    len(self.marker_set.leafes),
                    len(self.marker_set.covered),
                    len(self.marker_set.all_places),
                )
            )

    def LCA(self, leafes, debug=False):
        # reduce leafes to possible leafes
        logging.debug("Computing LCA from placement")
        leafes = list(set(leafes) & set(self.t.get_leaf_names()))

        if debug and len(leafes) == 1:
            # there is a bug with one leaf
            return self.t.get_tree_root().search_nodes(name=leafes[0])[0].up
        return self.t.get_common_ancestor(leafes)

    def node_taxid(self, node, info, return_lng=False):
        leafes = node.get_leaf_names()
        # reduce to possible keys, removing placements
        leafes = set(leafes) & set(list(info.keys()))
        lngs = [info[acc] for acc in leafes]

        # determine LCA
        if len(lngs) > 0:
            # determine LCA of all
            common = set(lngs[0])
            for lng in lngs[1:]:
                common = common & set(lng)
            shared = [x for x in lngs[0] if x in common]
            taxid = shared[-1]
            logging.debug("LCA taxid determined as {}".format(taxid))
            if return_lng:
                return shared
            else:
                return taxid

        else:
            # return 1 as default
            return 1

    def _find_best_ncbi_set(
        self,
        places,
        scmg,
        taxinfo,
        set_species=3,
        min_set_size=10,
        min_prevalence=95,
        set_atmost=500,
        sort_using="prevalence",
        training=False,
    ):

        # first determine LCA
        logging.debug("Using taxid to select marker gene sets")
        # loading taxid
        info = load_tax_info(taxinfo)
        # get lac_lineage
        LCA_lng = self.node_taxid(self.LCA(places), info, return_lng=True)

        all_sets = []
        mp = 100
        while mp >= min_prevalence:
            logging.debug("Searching marker gene sets at {}% prevalence".format(mp))
            # determine set at this position
            for taxid in LCA_lng[::-1]:
                # for each taxid try to define a set
                # first identify all possible leafes
                accs = [acc for acc, lng in info.items() if taxid in lng]
                accs = set(accs) & set(self.t.get_leaf_names())

                # obtain a marker gene set
                scmg_set = self._calc_set(
                    self.t.get_tree_root(),
                    places,
                    scmg,
                    set_species,
                    min_set_size,
                    mp,
                    set_atmost,
                    leafes=accs,
                )

                if scmg_set is not None:
                    # add gerneral info
                    # scmg_set.all_leafes = all_leafes
                    scmg_set.all_places = places
                    all_sets.append(scmg_set)

            mp = mp - 0.5
        return self._select_set(all_sets, sort_using, training)

    def _find_best_set(
        self,
        places,
        scmg,
        set_species=3,
        min_set_size=10,
        min_prevalence=95,
        set_atmost=500,
        training=False,
        sort_using="prevalence",
    ):

        all_sets = []
        # starting set search at root
        node = self.t.get_tree_root()
        # leafes = set(node.get_leaf_names()) - set(places)
        leafes = self.known_leafes & set(node.get_leaf_names())
        if logging.DEBUG >= logging.root.level:
            placements = set(node.get_leaf_names()) - leafes
            print(self.LCA(placements, debug=True))

        # all_leafes = leafes.copy()
        # move from root to leaf
        logging.debug("searching for all possible marker gene sets")
        # while len(leafes) >= set_species:
        # set marker gene prevalence to 100%
        mp = 100
        while mp >= min_prevalence:
            logging.debug("Searching marker gene sets at {}% prevalence".format(mp))
            # determine set at this position
            for node in self.t.traverse():
                # skip sets with not enough species diversity
                # leafes = set(node.get_leaf_names()) - set(places)
                leafes = self.known_leafes & set(node.get_leaf_names())
                if len(leafes) < set_species and training is False:
                    continue

                # obtain a marker gene set
                scmg_set = self._calc_set(
                    node,
                    places,
                    scmg,
                    set_species,
                    min_set_size,
                    mp,
                    set_atmost,
                    training=training,
                )

                if scmg_set is not None:
                    # add gerneral info
                    # scmg_set.all_leafes = all_leafes
                    scmg_set.all_places = places
                    all_sets.append(scmg_set)

            # decrease prevalence by fixed amount
            mp = mp - 0.5

        return self._select_set(all_sets, sort_using, training)

    def _select_set(self, all_sets, sort_using="lm", training=False):
        # check that we have at least 1 useable set
        if len(all_sets) == 0:
            logging.error("Could not identify a single suitable marker gene set")
            return None
        # now we can rank the sets and return the best one
        # sort_using = "lprevalenec"
        if sort_using == "lm" or sort_using.startswith("l"):
            logging.debug("Using linear model score to choose best set.")
            sorted_sets = sorted(all_sets, key=lambda x: x.score, reverse=False)
        elif sort_using == "prevalence" or sort_using.startswith("p"):
            logging.debug("Using marker gene prevalence to choose best set.")
            sorted_sets = sorted(
                all_sets,
                key=lambda x: (x.covered, x.prevalence, x.nspecies, x.nprofiles),
                reverse=True,
            )
        elif sort_using == "species" or sort_using.startswith("s"):
            logging.debug("Using marker gene species spread to choose best set.")
            sorted_sets = sorted(
                all_sets,
                key=lambda x: (x.covered, x.nspecies, x.prevalence, x.nprofiles),
                reverse=True,
            )
        elif sort_using == "best_guess" or sort_using.startswith("b"):
            logging.debug("Using marker gene species spread to choose best set.")
            # reduce to sets that cover at least 80% of the most placed placements
            most_placed = max([len(x.covered) for x in all_sets])
            red_sets = [x for x in all_sets if len(x.covered) >= 0.8 * most_placed]
            sorted_sets = sorted(red_sets, key=lambda x: (x.best_guess), reverse=True)
        else:
            logging.error("No valid sorting algorythm specified, results are unsorted!")
            raise NotImplementedError("No valid sorting method")

        if training:
            return sorted_sets
        else:
            return sorted_sets[0]

    # @profile
    def _calc_set(
        self,
        node,
        places,
        scmg,
        set_species=5,
        set_size=50,
        set_prevalence=100,
        set_atmost=500,
        training=False,
        leafes=None,
    ):
        """
        given a node, report number of leafes and the set at this location
        """
        if leafes is None and node is not None:
            # leafes = set(node.get_leaf_names()) - set(places)
            leafes = self.known_leafes & set(node.get_leaf_names())
        sets = []
        for leaf in leafes:
            try:
                sets.append(scmg[leaf])
            except KeyError:
                logging.warning(
                    "Database missing markes for '{}'. This should not be the case. Make sure the database is not corrupted".format(
                        leaf
                    )
                )
        placements = set(node.get_leaf_names()) & set(places)
        n_placements = len(placements)
        # quickly terminate if no placements are in this potential location
        if n_placements == 0:
            return None
        s = percentage_sets(sets, set_prevalence, set_atmost)

        if logging.DEBUG >= logging.root.level:
            logging.debug(
                "Looking at set with size {} and {} leafes covering {}".format(
                    len(s), len(leafes), n_placements
                )
            )

        if (len(s) >= set_size and n_placements > 0) or training is True:
            logging.debug(
                "Set size {} with {} leafes, covering {}/{} placements".format(
                    len(s), len(leafes), n_placements, len(places)
                )
            )
            return markerset(leafes, placements, s, set_prevalence)
        else:
            return None
