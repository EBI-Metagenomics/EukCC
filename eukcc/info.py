import os
import json
from eukcc import base
from eukcc.base import log


class eukinfo():
    def __init__(self, dirname):
        self.dirname = dirname
        v = self.checkForFiles(dirname)
        
        self.placementHMMs = os.path.join(self.dirname, "hmms/concat.hmm")

    def checkForFiles(self, dirname):
        required = ["profile.list", "refpkg", "hmms/concat.hmm"]
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
        
    
    
    

