import os
from eukcc.info import eukinfo
from eukcc.base import log
from eukcc.exec import checkDependencies
from eukcc.exec import hmmer
from eukcc.exec import hmmpress
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
    def __init__(self, fastapath, configdir, 
                 outfile = "eukcc.tsv", 
                 threads = None,
                 outdir = None, 
                 place = None, verbose = True, force = None, 
                 isprotein = None, bedfile = None):
        # check config dir
        self.config = eukinfo(configdir)
        self.cfg = self.config.cfg
        # update config with function params
        self.cfg = updateConf(self.cfg, "outdir", outdir)
        self.cfg = updateConf(self.cfg, "force", force)
        self.cfg = updateConf(self.cfg, "threads", threads)
        self.cfg = updateConf(self.cfg, "verbose", verbose)
        self.cfg = updateConf(self.cfg, "isprotein", isprotein)
        self.cfg = updateConf(self.cfg, "place", place)
        self.cfg = updateConf(self.cfg, "outfile", outfile)

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
            self.placed = self.place(proteinfaa, bedfile)
        else:
            print("check if placement can be found in tree")
            print("if so, use that and run next step")
        
        # compute completeness and contamination
        hmmfile = self.concatHMM(self.placed)
        hits = self.runPlacedHMM(hmmfile, proteinfaa, bedfile)
        
        outputfile = os.path.join(self.cfg['outdir'], self.cfg['outfile'])
        self.estimate(hits, outputfile, self.placed)
        
    def checkIO(self, fastapath, outdir):
        # create outdir if not exists
        file.isdir(self.cfg['outdir'])

        print("IO check not yet implemented")
        return(False)
    
    def estimate(self, hits, outfile, placements):
        hit = {}
        r = base.readTSV(hits)
        # count profile hits
        for row in r:
            if row['profile'] not in hit.keys():
                hit[row['profile']] = 0
            else:
                hit[row['profile']] += 1
        doubles = set([n for k,v  in hit.items() if v > 1])
        # now we can estimate completeness and contamination for each placement
        for i in range(len(placements)):
            s = self.readSet(placements[i]['path'])
            # completeness is the overap of both sets 
            cmpl = len(s & set(hit.keys()))/len(s)
            cont = len(doubles & s)/len(s)
            placements[i]['completeness'] = round(cmpl * 100, 2)
            placements[i]['contamination'] = round(cont * 100, 4)
        
        print(outfile)
        k = ["n", "ngenomes", "completeness", "contamination", "tax_id", "cover"]
        with open(outfile, "w") as f:
            f.write("{}\n".format("\t".join(k)))
            for p in placements:
                f.write("\t".join([str(p[key]) for key in k]+["\n"]))
                    
            
        return(placements)
    
    def readSet(self, p):
        profiles = []
        with open(p) as f:
            for line in f:
                profiles.append(line.strip())
        return(set(profiles))
                    
    
    def concatHMM(self, places):
        profiles = []
        for p in places:
            with open(p['path']) as f:
                for line in f:
                    profiles.append(line.strip())
        # create all paths for all hmms
        hmmerpaths = [os.path.join(self.cfg['PANTHER'], "books", profile, "hmmer.hmm") for profile in profiles]
        
        # create a dir for this
        hmmdir = os.path.join(self.cfg['outdir'],"workfiles","hmmer", "estimations")
        file.isdir(hmmdir)
        hmmconcat = os.path.join(hmmdir, "all.hmm")
        # concatenate
        log("{} hmm profiles need to be used for estimations".format(len(profiles)), self.cfg['verbose'])
        log("Concatenating hmms, this might take a while (IO limited)", self.cfg['verbose'])
        #print("Unkomment this here for production!")
        hmmconcat = base.concatenate(hmmconcat, hmmerpaths)
        # press
        log("Pressing hmms", self.cfg['verbose'])
        hp = hmmpress("hmmpress", hmmconcat, None)
        hp.run()
        return(hmmconcat)

    def runPlacedHMM(self, hmmfile, proteinfaa, bedfile):
        # run hmmer and strip down

        # define output files
        hmmDir = os.path.join(self.cfg['outdir'],"workfiles","hmmer", "estimations")
        file.isdir(hmmDir)
        hmmOut = os.path.join(hmmDir, "placement.tsv")
        hmmOus = os.path.join(hmmDir, "placement.out")
        hitOut = os.path.join(hmmDir, "hits.tsv")
        
        h = hmmer("hmmsearch", proteinfaa, hmmOut)
        if h.doIneedTorun(self.cfg['force']):
            log("Running hmmer for chosen locations", self.cfg['verbose'])
            h.run(hmmOus, hmmfiles = hmmfile, 
                  cores = self.cfg['threads'])
            # clean hmmer outpout
            log("Processing Hmmer results", self.cfg['verbose'])
            hitOut = h.clean(hmmOut, bedfile, hitOut, self.cfg['mindist'])
        return(hitOut)
        
    
    def gmes(self, fasta):
        """
        predict proteins using gmes
        """
        gmesDir = os.path.join(self.cfg['outdir'],"workfiles","gmes")
        file.isdir(gmesDir)
        gmesOut = os.path.join(gmesDir, "proteins.faa")
        gtffile = os.path.join(gmesDir, "genemark.gtf")

        g = gmes("runGMES", fasta, gmesOut)
        if g.doIneedTorun(self.cfg['force']):
            g.run(cores = self.cfg['threads'])
        
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
            h.run(hmmOus, hmmfiles = self.config.placementHMMs,
                  cores = self.cfg['threads'])
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
        pplaceOutReduced = os.path.join(placerDir, "placementReduced.jplace")
        file.isdir(placerDirTmp)
        
        # pplacer
        pp = pplacer("pplacer",pplaceAlinment, pplaceOut)
        if pp.doIneedTorun(self.cfg['force']):
            log("Running pplacer", self.cfg['verbose'])
            pp.prepareAlignment(hitOut, os.path.join(self.config.dirname, "profile.list"), fasta, 
                                self.config, self.cfg,  placerDirTmp )
            pp.run(os.path.join(self.config.dirname, "refpkg", "concat.refpkg"),
                   cores = self.cfg['threads'])
        
        # reduce placements to the placements with at least posterior of p
        pplaceOutReduced = pp.reduceJplace(pplaceOut, pplaceOutReduced, self.cfg['minPlacementLikelyhood'])
        
        # run TOG to get a tree
        togTree = os.path.join(placerDir, "placement.tree")
        tg = tog("guppy", pplaceOutReduced, togTree)
        if tg.doIneedTorun(self.cfg['force']):
            log("Fetching pplacer tree", self.cfg['verbose'])
            tg.run()
        
        log("Getting best placement", self.cfg['verbose'])
        # now we can place the bin using the tree
        t = treelineage.treeHandler(togTree)
        t2 = treelineage.treeHandler(self.config.tree)
        orignialleaves = t2.leaves()
        sets = self.getSets()
        placements = t.getPlacement(self.cfg['placementMethod'], sets, orignialleaves, self.cfg['nPlacements'], self.cfg['minSupport'])
        log("Done placeing", self.cfg['verbose'])
        return(placements)
        
    
    def getSets(self):
        setinfo = os.path.join(self.config.dirname, "sets", "setinfo.csv")
        # load sets and reduce to sets matching parameters 
        sets = []
        with open(setinfo) as f:
            cols = []
            for line in f:
                line = line.strip()
                l = line.split(",")
                if len(cols) == 0:
                    cols = l
                    continue
                n = {}
                for k, v in zip(cols, l):
                    n[k] = v
                if int(n["ngenomes"]) >= self.cfg['minGenomes'] and int(n['n']) >= self.cfg['minProfiles']:
                    sets.append(n)
                    
        return(sets)
                
        
        
        
