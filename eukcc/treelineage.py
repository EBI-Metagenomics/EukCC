#!/usr/bin/env python
import operator
from ete3 import Tree
from ete3 import NCBITaxa
ncbi = NCBITaxa()


class treeHandler():
    def __init__(self, tree = None, root = False, annotate = True):
        self.loadTree(tree)
        if root:
            self.root()
        if annotate:
            self.annotateTree()
        
    def loadTree(self, tree):
        try:
            self.t = Tree(tree)
        except:
            self.t = Tree(tree, format = 1)
        
    def annotateTree(self):
        i = 0
        for node in self.t.get_tree_root().traverse():
            if node.name == "":
                node.name = "node{}".format(i)
                i += 1
    
    def get_lineage(self, leave, istaxid = True):
        try:
            if istaxid:
                n = self.t.search_nodes(taxid=int(leave))[0]
            else:
                n = self.t.search_nodes(name=leave)[0]
        except:
            print("Node {} not found".format(leave))
            return(False)
        l = []
        while True:
            try:
                l.append(n.name)
                n = n.up
            except:
                break
        l.reverse()
        return(l)
    
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
        return(self.t.get_tree_root().get_leaf_names())
    
    def children(self, nodename ):
        '''
        find and rturn all leafs belpongin as chidlren
        to a node
        '''
        nn = self.t.search_nodes(name=nodename)[0].get_leaf_names()
        #for node in self.t.search_nodes(name=nodename)[0].iter_descendants():
        #    if node.is_leaf():
        #        nn.append(node.name)
        return(nn)

        
    def write(self, file):
        self.t.write(format=1, outfile=file)
    
    def getPlacement(self, mode, sets, nonplacements, nplacements = 2, atleast = 1):
        """
        function to find the set that has the most support given either LCA (default)
        or HPA placement .
        returns list of dict
        """
        placements = set(self.t.get_tree_root().get_leaf_names()) - set(nonplacements)
        remaining = placements
        results = []
        while nplacements > 0:
            for i in range(len(sets)):
                covering = set(self.children(sets[i]['tax_id'])) & remaining
                sets[i]['cover'] = len(covering)
                sets[i]['covering'] = covering
                sets[i]['nPlacements'] = len(placements)
            
            # get the one with the best coverage sorted by what we need now
            # so we need to sort the list of dicts with different keys
            if mode == "HPA":
                sets.sort(key=operator.itemgetter("n"), reverse=True)
                sets.sort(key=operator.itemgetter("ngenomes"), reverse=True)
            else:
                sets.sort(key=operator.itemgetter("ngenomes"), reverse=True)
                sets.sort(key=operator.itemgetter("n"), reverse=True)
            sets.sort(key=operator.itemgetter("cover"), reverse=True)
            
        
            # only retain if at least N placements
            if sets[0]['cover'] >= atleast:
                results.append(sets[0].copy())
                # remove remaining
                remaining = remaining - sets[0]['covering']
            
            nplacements -= 1
            
        return(results)
        
    
    def getLCA(self, placements, nodes):
        return



if __name__ == "__main__":
    n = treeHandler("/hps/nobackup2/production/metagenomics/saary/markergenes2/trees/trees/mafft/fasttree/concat_midpointrooted.tree")
    print(n.get_lineage("GCA_001625265.1", False))
    n.write("/hps/nobackup2/production/metagenomics/saary/markergenes2/trees/trees/mafft/fasttree/concat_rooted.tree")


# In[ ]:





# In[126]:





# In[10]:





# In[ ]:




