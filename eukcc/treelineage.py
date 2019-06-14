#!/usr/bin/env python
# coding: utf-8

# In[14]:


from ete3 import Tree
from ete3 import NCBITaxa
ncbi = NCBITaxa()


# In[113]:


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
        
    def children(self, nodename ):
        '''
        find and rturn all leafs belpongin as chidlren
        to a node
        '''
        nn = []
        for node in self.t.search_nodes(name=nodename)[0].iter_descendants():
            if node.is_leaf():
                nn.append(node.name)
        return(nn)

        
    def write(self, file):
        self.t.write(format=1, outfile=file)


# In[116]:


if __name__ == "__main__":
    n = treeHandler("/hps/nobackup2/production/metagenomics/saary/markergenes2/trees/trees/mafft/fasttree/concat_midpointrooted.tree")
    print(n.get_lineage("GCA_001625265.1", False))
    n.write("/hps/nobackup2/production/metagenomics/saary/markergenes2/trees/trees/mafft/fasttree/concat_rooted.tree")


# In[ ]:





# In[126]:





# In[10]:





# In[ ]:




