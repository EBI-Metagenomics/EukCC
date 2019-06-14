import os
from eukcc.info import eukinfo
from eukcc.base import log
from eukcc.exec import checkDependencies
from eukcc.exec import hmmer
from eukcc.exec import gmes
from eukcc.exec import pplacer
from eukcc.exec import tog
from eukcc import base
from eukcc import treelineage

from eukcc.fileoperations import file

dep = {
       "hmmer": "hmmsearch",
       "runGMES": "runGMES",
       "GeneMarkES": "gmes_petap.pl"
      }


def updateConf(cfg, k, v):
    if v is not None:
        cfg[k] = v
    return(cfg)
    

class eukcc():
    def __init__(self, fastapath, configdir, outdir = None, place = None, verbose = True, force = None, isprotein = None, bedfile = None):
        # check config dir
        self.config = eukinfo(configdir)
        self.cfg = self.config.cfg
        # update config with function params
        self.cfg = updateConf(self.cfg, "outdir", outdir)
        self.cfg = updateConf(self.cfg, "force", force)
        self.cfg = updateConf(self.cfg, "verbose", verbose)
        self.cfg = updateConf(self.cfg, "isprotein", isprotein)
        self.cfg = updateConf(self.cfg, "place", place)

        # check if we can read and write
        self.checkIO(fastapath, outdir)
        # check for dependencies
        checkDependencies(dep)
        
        # skip gene predition if this is already protein sequences
        if not isprotein:
            # run gmes
            proteinfaa, bedfile = self.gmes(fastapath)
        else:
            proteinfaa = fastapath
        
        # place
        if place is None:
            self.place(proteinfaa, bedfile)
        else:
            print("check if placement can be found in tree")
            print("if so, use that and run next step")
        
        # compute completeness and contamination
        
    def checkIO(self, fastapath, outdir):
        # create outdir if not exists
        file.isdir(self.cfg['outdir'])

        print("IO check not yet implemented")
        return(False)
    
    
    def gmes(self, fasta):
        """
        predict proteins using gmes
        """
        gmesDir = os.path.join(self.cfg['outdir'],"workfiles","gmes")
        file.isdir(gmesDir)
        gmesOut = os.path.join(gmesDir, "proteins.faa")
        gtffile = gmesOut[:-4] + ".gtf"

        g = gmes("runGMES", fasta, gmesOut)
        if g.doIneedTorun(self.cfg['force']):
            g.run()
        
        # make a bed file from GTF
        bedf = os.path.join(gmesDir, "proteins.bed")
        if self.cfg['force'] or file.isnewer(gtffile, bedf):
            
            bedf = base.gmesBED(gtffile, bedf)
        print(gtffile)
        return(gmesOut, bedf)
        
        
    
    def place(self, fasta, bedfile):
        """
        main function to place a bin in the tree.
        will subsequently run hmmer 
        """
        # test if we can open the input files first
        if not base.exists(fasta):
            print("Could not open fasta file")
            return(False)
        if not base.exists(bedfile):
            print("Could not open bedfile file")
            return(False)
        
        
        # define output files
        hmmDir = os.path.join(self.cfg['outdir'],"workfiles","hmmer")
        file.isdir(hmmDir)
        hmmOut = os.path.join(hmmDir, "placement.tsv")
        hmmOus = os.path.join(hmmDir, "placement.out")
        hitOut = os.path.join(hmmDir, "hits.tsv")
        
        # run hmmer if forced or input newer than output            
        h = hmmer("hmmsearch", fasta, hmmOut)
        if h.doIneedTorun(self.cfg['force']):
            log("Gonna place the bin", self.cfg['verbose'])
            h.run(hmmOus, hmmfiles = self.config.placementHMMs)
            # clean hmmer outpout
            log("Processing Hmmer results", self.cfg['verbose'])
            hitOut = h.clean(hmmOut, bedfile, hitOut, self.cfg['mindist'])
        else:
            log("Bin was already placed", self.cfg['verbose'])
        
        # pplacer paths
        placerDir = os.path.join(self.cfg['outdir'],"workfiles","pplacer")
        placerDirTmp = os.path.join(placerDir, "tmp")
        pplaceAlinment = os.path.join(placerDir, "horizontalAlignment.fasta")
        pplaceOut = os.path.join(placerDir, "placement.jplace")
        file.isdir(placerDirTmp)
        
        # pplacer
        pp = pplacer("pplacer",pplaceAlinment, pplaceOut)
        if pp.doIneedTorun(self.cfg['force']):
            log("Running pplacer", self.cfg['verbose'])
            pp.prepareAlignment(hitOut, os.path.join(self.config.dirname, "profile.list"), fasta, 
                                self.config, self.cfg,  placerDirTmp )
            pp.run(os.path.join(self.config.dirname, "refpkg", "concat.refpkg"))
        
        # run TOG to get a tree
        togTree = os.path.join(placerDir, "placement.tree")
        tg = tog("guppy", pplaceOut, togTree)
        if tg.doIneedTorun(self.cfg['force']):
            log("Fetching pplacer tree", self.cfg['verbose'])
            tg.run()
            
        # now we can place the bin using the tree
        t = treelineage.treeHandler(togTree)
        
        
        
        
