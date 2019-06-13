import os
from eukcc.info import eukinfo
from eukcc.base import log
from eukcc.exec import checkDependencies
from eukcc.exec import hmmer
from eukcc.exec import gmes

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
    def __init__(self, fastapath, configdir, outdir = None, place = None, verbose = True, force = None, isprotein = None):
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
            proteinfaa = self.gmes(fastapath)
        else:
            proteinfaa = fastapath
        
        # place
        if place is None:
            self.place(proteinfaa)
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

        g = gmes("runGMES", fasta, gmesOut)
        if g.doIneedTorun(self.cfg['force']):
            g.run()
            
        return(gmesOut)
        
        
    
    def place(self, fasta):
        """
        main function to place a bin in the tree.
        will subsequently run hmmer 
        """
        
        # define output files
        
        hmmDir = os.path.join(self.cfg['outdir'],"workfiles","hmmer")
        file.isdir(hmmDir)
        hmmOut = os.path.join(hmmDir, "placement.tsv")
        hmmOus = os.path.join(hmmDir, "placement.out")
        
        # run hmmer if forced or input newer than output            
        h = hmmer("hmmsearch", fasta, hmmOut)
        if h.doIneedTorun(self.cfg['force']):
            log("Gonna place the bin", self.cfg['verbose'])
            h.run(hmmOus, hmmfiles = self.config.placementHMMs)
            # clean hmmer outpout
        else:
            log("Bin was already placed", self.cfg['verbose'])

        
        
        
        
        
