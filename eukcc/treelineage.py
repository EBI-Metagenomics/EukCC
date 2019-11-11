#!/usr/bin/env python
import operator
import os
import json
import logging
from ete3 import Tree, NodeStyle, TreeStyle
from ete3 import NCBITaxa

ncbi = NCBITaxa()


def RGB_to_hex(RGB):
    """ [255,255,255] -> "#FFFFFF" """
    # Components need to be integers for hex to make sense
    RGB = [int(x) for x in RGB]
    return "#" + "".join(["0{0:x}".format(v) if v < 16 else "{0:x}".format(v) for v in RGB])


class treeHandler:
    def __init__(self, tree=None, root=False, annotate=True):
        self.loadTree(tree)
        if root:
            self.root()
        if annotate:
            self.annotateTree()
        else:
            # always name root node0
            self.t.get_tree_root().name = "node0"

    def loadTree(self, tree):
        try:
            self.t = Tree(tree)
        except:
            self.t = Tree(tree, format=1)

    def annotateTree(self):
        i = 0
        for node in self.t.get_tree_root().traverse():
            if node.name == "":
                node.name = "node{}".format(i)
                i += 1

    def get_lineage(self, leave, istaxid=True):
        try:
            if istaxid:
                n = self.t.search_nodes(taxid=int(leave))[0]
            else:
                n = self.t.search_nodes(name=leave)[0]
        except:
            print("Node {} not found".format(leave))
            return False
        l = []
        while True:
            try:
                l.append(n.name)
                n = n.up
            except:
                break
        l.reverse()
        return l

    def addTaxids(self, GCAs, taxids):
        for GCA, tax in zip(GCAs, taxids):
            try:
                n = self.t.search_nodes(name=GCA)[0]
                n.add_feature("taxid", int(tax))
            except:
                continue

    def root(self):
        # midpoitn rooting
        R = self.t.get_midpoint_outgroup()
        # and set it as tree outgroup
        self.t.set_outgroup(R)

    def leaves(self):
        return self.t.get_tree_root().get_leaf_names()

    def children(self, nodename):
        """
        find and rturn all leafs belpongin as chidlren
        to a node
        """
        try:
            nn = self.t.search_nodes(name=nodename)[0].get_leaf_names()
        except:
            print(f"Could not find any children for node {nodename}")
            nn = []
        return nn

    def write(self, file):
        self.t.write(format=1, outfile=file)

    def getPlacement(self, mode, sets, originaltree, nplacements=2, atleast=1, maximum=3, debug=False):
        """
        function to find the set that has the most support given either LCA (default)
        or HPA placement .
        returns list of dict
        """
        nonplacements = originaltree.leaves()
        # get all placements by subtsracting known leaves from all leave names
        placements = set(self.t.get_tree_root().get_leaf_names()) - set(nonplacements)
        remaining = placements
        results = []
        # while we have placements to cover, search for sets
        while nplacements > 0:
            for i in range(len(sets)):
                covering = set(self.children(sets[i]["node"])) & remaining
                sets[i]["cover"] = len(covering)
                sets[i]["covering"] = covering
                sets[i]["nPlacements"] = len(placements)
            # get the one with the best coverage sorted by what we need now
            # so we need to sort the list of dicts with different keys
            if mode == "HPA":
                # sets.sort(key=operator.itemgetter("n"), reverse=True)
                sets.sort(key=operator.itemgetter("ngenomes"), reverse=True)
            elif mode == "LCA":
                # sets.sort(key=operator.itemgetter("n"), reverse=True)
                sets.sort(key=operator.itemgetter("ngenomes"), reverse=False)
            else:
                # we dont know what to do
                logging.warning("Mode not know", mode)
                exit()
            # sort now by cover, keeping the underlying order of genomes in case
            # several sets cover the same amount of profiles
            sets.sort(key=operator.itemgetter("cover"), reverse=True)

            # only retain if at least N placements
            i = 0
            # in debug mode we want to return all best placements, not just the LCA or HPA
            while i < maximum and sets[i]["cover"] >= atleast:
                # break if new set has less more than 1 less covers than the best
                if i > 0 and (results[0]["cover"] - sets[i]["cover"]) > 1:
                    break

                # save set as a result
                results.append(sets[i].copy())
                # add neighbours
                sisters = []
                # last result item we track which placmeent where covered
                for place in results[-1]["covering"]:
                    # we search the sisters of the placement loaction
                    si = self.t.search_nodes(name=place)[0].get_sisters()
                    for sis in si:
                        # and fetch the names of all leaves
                        sisters.extend(sis.get_leaf_names())
                # reducing this to only keep GCA numbers and not pplacer leaves
                sisters = set(sisters) - placements
                # save in list
                results[-1]["sisters"] = sisters
                # remove remaining
                remaining = remaining - sets[i]["covering"]

                # count up i, usually
                i = i + 1
                logging.debug(results[-1])

            nplacements -= 1

        return results

    def loadInfo(self, togjson):
        info = {}
        with open(togjson) as json_file:
            j = json.load(json_file)
            fields = j["fields"]
            # loop placements
            for placement in j["placements"]:
                nm = placement["nm"][0][0]
                ps = placement["p"]
                kp = []
                for p in ps:
                    d = {k: v for k, v in zip(fields, p)}

                info[nm] = d
        return info

    def plot(self, placement, togjson, outdir, cfg):
        """
        plot a plcement in the tree
        show all pplacer placements and the LCA and HCA node 
        as well as the inferred lineage
        """
        # with no X display this needs to be set
        os.environ["QT_QPA_PLATFORM"] = "offscreen"
        info = self.loadInfo(togjson)

        no = 0
        for LCAp, HPAp in zip(placement["LCA"], placement["HPA"]):

            plotpath = os.path.join(outdir, f"tree_{no}.png")

            # make shallow copy
            t = self.t

            LCA = LCAp["node"]
            HPA = HPAp["node"]
            # define basic tree style
            ts = TreeStyle()
            # hide leave names
            ts.show_leaf_name = False

            # circular tree
            ts.mode = "c"
            ts.rotation = 210
            ts.arc_start = 0  # 0 degrees = 3 o'clock
            ts.arc_span = 350

            highlightsize = 50
            nodesize = 10
            # define styles for special nodes
            # at the moment hard coded, but could be accesible for the user
            # PTHR style
            nstyle = NodeStyle()
            nstyle["fgcolor"] = "red"
            nstyle["size"] = highlightsize

            # LCA style
            LCAstyle = NodeStyle()
            LCAstyle["fgcolor"] = "green"
            LCAstyle["bgcolor"] = "DarkSeaGreen"
            LCAstyle["size"] = highlightsize

            # HPA style
            HPAstyle = NodeStyle()
            HPAstyle["fgcolor"] = "blue"
            HPAstyle["bgcolor"] = "LightSteelBlue"
            HPAstyle["size"] = highlightsize

            # default node
            defaultStyle = NodeStyle()
            defaultStyle["fgcolor"] = "gray"
            defaultStyle["size"] = nodesize

            for n in t.traverse():
                if n.name.startswith("PTHR"):
                    # set color based on posterior prob:
                    x = (info[n.name]["post_prob"] - cfg["minPlacementLikelyhood"]) / (
                        1 - cfg["minPlacementLikelyhood"]
                    )
                    c = [int(x * 255), (1 - x) * 120, 0]
                    he = RGB_to_hex(c)
                    # define back color of locations
                    nstyle["bgcolor"] = he
                    n.set_style(nstyle)

                elif n.name == LCA:
                    n.set_style(LCAstyle)
                elif n.name == HPA:
                    n.set_style(HPAstyle)
                else:
                    n.set_style(defaultStyle)

            # plot to disk
            _ = t.render(plotpath, w=320, units="mm", tree_style=ts)
            no = no + 1
