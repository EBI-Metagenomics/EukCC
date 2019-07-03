import os
import json
from eukcc import base
from eukcc.base import log
import yaml


defaults = {"verbose": True,
            "outfile": "eukcc.tsv",
            "isprotein": False,
            "outdir": ".",
            "place": None,
            "hmm": False,
            "force": False,
            "threads": 1,
            "placementMethod": "LCA",
            "minProfiles": 20,
            "minGenomes": 3,
            "minSupport": 2,
            "nPlacements": 1,
            "noplace": False,
            "cleanfasta": True,
            "PANTHER": "/homes/saary/data/databases/panther/PANTHER_14.1/PANTHER14.1",
            "minPlacementLikelyhood": 0.5,
            "evalue": 1e-5,
            "mindist": 2000}


class eukinfo():
    def __init__(self, dirname):
        self.dirname = dirname
        v = self.checkForFiles(dirname)
        # define location of placement HMMs
        self.placementHMMs = os.path.join(self.dirname, "hmms/concat.hmm")
        self.tree = self.pkgfile("concat.refpkg", "tree")
        
        # define deaults and load config if any
        self.cfg = defaults
        self.loadConfig()

    def checkForFiles(self, dirname):
        required = ["profile.list", "refpkg", "hmms/concat.hmm", "sets/setinfo.csv"]
        for f in required:
            p = os.path.join(dirname, f)
            if not base.exists(p):
                print("Configuartion folder does not contain: {}".format(f))
                return(False)
        return(True)
    
    def loadConfig(self):
        """
        See if config.yaml can be found and if so
        load it and iverwrite defaults
        """
        cp = os.path.join(self.dirname, "config.yaml")
        if not base.exists(cp):
            return
        with open(cp) as f:
            cfg = yaml.load(f)
            for k,v in cfg.items():
                self.cfg[k] = v
    
    
    def pkgfile(self, name, t):
        """
        get a file path for a refpkg package
        """
        info = self.readInfo(name)
        p = os.path.join(self.dirname, "refpkg", name, info['files'][t])
        if base.exists(p):
            return(p)
        else:
            log("Could not find: {}".format(p))
        
    def readInfo(self, name):
        p = os.path.join(self.dirname, "refpkg", name, "CONTENTS.json")
        # raise error if we cant find the file
        if not base.exists(p):
            log("Could not find {}".format(p))
            return
        # read and return json
        with open(p) as json_file:  
            j = json.load(json_file)
            return(j)
        
    
    
    

