import os
import json
import logging
from eukcc import base
from eukcc.base import log


old_defaults = {"verbose": True,
            "debug": False,
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
            "nPlacements": 2,
            "noplace": False,
            "lineage": "limited", # limited or full
            "cleanfasta": True,
            "minPlacementLikelyhood": 0.8,
            "evalue": 1e-5,
            "trainingEvalue": 10,
            "training": False,
            "mindist": 2000,
            "steps": {"gmes": False,
                      "findprots": False,
                      "pplacer": False,
                      "runHmmer": False,
                      "estimatedCompleteness": False}}

# default magic numbers are saved here
# please keep it sorted alphabetically and coument
defaults = {
            "evalue": 1e-5,        # not used evalue
            "nEvals": 3,           # per placement show top n inferrals
            "minGenomes": 3,       # minimal number of genomes in set
            "minProfiles": 20,     # minimal profiles in set
            "minSupport": 2,       # minimal profiles to support set
            "trainingEvalue": 10   # evalue used in training mode
            }


class eukinfo():
    def __init__(self, options):
        self.setConfig(options)


    def setConfig(self, options):
        self.cfg = defaults
        # loop over options and import into our configs
        for k, v in vars(options).items():
            self.cfg[k] = v
        
        # placements Methis
        if options.HPA:
            self.cfg['placementMethod'] = "HPA"
        else:
            self.cfg['placementMethod'] = "LCA"
        # define location of placement HMMs
        self.placementHMMs = os.path.join(options.db, "hmms/concat.hmm")
        self.tree = self.pkgfile("concat.refpkg", "tree")
        

    def checkForFiles(self, dirname):
        required = ["profile.list", 
                    "refpkg", 
                    "hmms/concat.hmm",
                    "sets/setinfo.csv"]
        for f in required:
            p = os.path.join(dirname, f)
            if not base.exists(p):
                print("Configuartion folder does not contain: {}".format(f))
                return(False)
        return(True)

    def pkgfile(self, name, t):
        """
        get a file path for a refpkg package
        """
        info = self.readInfo(name)
        p = os.path.join(self.cfg['db'], "refpkg", name, info['files'][t])
        if base.exists(p):
            return(p)
        else:
            log("Could not find: {}".format(p))
            exit()

    def readInfo(self, name):
        p = os.path.join(self.cfg['db'], "refpkg", name, "CONTENTS.json")
        # raise error if we cant find the file
        if not base.exists(p):
            log("Could not find {}".format(p))
            exit()
        # read and return json
        with open(p) as json_file:
            j = json.load(json_file)
            return(j)
