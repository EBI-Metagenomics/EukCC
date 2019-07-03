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
       "GeneMarkES": "gmes_petap.pl",
       "pplacer": "pplacer"
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
                 place = None, 
                 verbose = True, 
                 force = None, 
                 fplace = None,
                 cleanfasta = None,
                 isprotein = None, 
                 bedfile = None,
                 hmm = None,
                 noplace = None
                ):
        # check config dir
        self.config = eukinfo(configdir)
        self.cfg = self.config.cfg
        # update config with function params
        self.cfg = updateConf(self.cfg, "cleanfasta", cleanfasta)
        self.cfg = updateConf(self.cfg, "force", force)
        self.cfg = updateConf(self.cfg, "fplace", fplace)
        self.cfg = updateConf(self.cfg, "isprotein", isprotein)
        self.cfg = updateConf(self.cfg, "hmm", hmm)
        self.cfg = updateConf(self.cfg, "noplace", noplace)
        self.cfg = updateConf(self.cfg, "outdir", outdir)
        self.cfg = updateConf(self.cfg, "outfile", outfile)
        self.cfg = updateConf(self.cfg, "place", place)
        self.cfg = updateConf(self.cfg, "threads", threads)
        self.cfg = updateConf(self.cfg, "verbose", verbose)
        
        self.stopped = {"stopped": False,
                        "reason": ""}

        # check if we can read and write
        self.checkIO(fastapath, outdir)
        # check for dependencies
        checkDependencies(dep)
        
        # skip gene predition if this is already protein sequences
        if not isprotein and not self.stopnow():
            # run gmes
            proteinfaa, bedfile = self.gmes(fastapath)
        else:
            proteinfaa = fastapath        
        
        # run hmm file if we are asked to
        # this is needed during for training 
        if self.cfg['hmm'] and not self.stopnow():
            _a = self.runPlacedHMM(self.cfg['hmm'], proteinfaa, bedfile)
            
        
        if self.cfg['noplace']:
            self.stop("Stopping because we were told to")
        
        # place using pplacer and hmmer
        if place is None and not self.stopnow():
            self.placed = self.place(proteinfaa, bedfile)
        else:
            print("check if placement can be found in tree")
            print("if so, use that and run next step")
        
        # concat hmms for hmmer
        if not self.stopnow():
            hmmfile = self.concatHMM(self.placed)
        # run Hmmer for sets of placement
        if not self.stopnow():
            hits = self.runPlacedHMM(hmmfile, proteinfaa, bedfile)
        # estimate completeness and contamiantion
        if not self.stopnow():            
            outputfile = os.path.join(self.cfg['outdir'], self.cfg['outfile'])
            self.estimate(hits, outputfile, self.placed)
    
    def stopnow(self):
        '''
        check if we need to stop bc something failed
        '''
        return(self.stopped['stopped'])
    
    def stop(self, reason):
        """
        set the stop info after a step failed
        """
        self.stopped = {"stopped": True,
                        "reason": reason}
        log("Stopping because: {}".format(reason))
        
    def checkIO(self, fastapath, outdir):
        # create outdir if not exists
        file.isdir(self.cfg['outdir'])
        # check if input and output can be accessed
        print("Warning: IO check not yet implemented")
        return(False)
    
    def estimate(self, hits, outfile, placements):
        hit = {}
        r = base.readTSV(hits)
        # count profile hits
        for row in r:
            if row['profile'] not in hit.keys():
                hit[row['profile']] = 1
            else:
                hit[row['profile']] += 1
        
        singletons = set(hit.keys())
        multitons = set([k for k,v  in hit.items() if v > 1])

        # now we can estimate completeness and contamination for each placement
        for i in range(len(placements)):
            s = self.readSet(placements[i]['path'])
            # completeness is the overap of both sets 
            cmpl = len(singletons & s)/len(s)
            cont = len(multitons & s)/len(s)
            # make to percentage and round to 2 positions
            placements[i]['completeness'] = round(cmpl * 100, 2)
            placements[i]['contamination'] = round(cont * 100, 2)
        
        log("finished estimating", self.cfg['verbose'])
        
        # write to output file
        k = ["completeness", "contamination", "tax_id", "n", "ngenomes", "cover", "nPlacements"]
        with open(outfile, "w") as f:
            f.write("{}\n".format("\t".join(k)))
            for p in placements:
                f.write("\t".join([str(p[key]) for key in k]+["\n"]))
        
        log("Wrote estimates to: {}".format(outfile), self.cfg['verbose'])
        
        # done
        return(True)
    
    def readSet(self, p):
        profiles = []
        localpath = os.path.join(self.config.dirname, "sets", os.path.basename(p))
        with open(localpath) as f:
            for line in f:
                profiles.append(line.strip())
        return(set(profiles))        
    
    def concatHMM(self, places):
        profiles = []
        for p in places:
            localpath = os.path.join(self.config.dirname, "sets", os.path.basename(p['path']))
            with open(localpath) as f:
                for line in f:
                    profiles.append(line.strip())
        # create all paths for all hmms
        hmmerpaths = [os.path.join(self.cfg['PANTHER'], "books", profile, "hmmer.hmm") 
                      for profile in profiles]
        
        # create a dir for this
        hmmdir = os.path.join(self.cfg['outdir'],"workfiles","hmmer", "estimations")
        file.isdir(hmmdir)
        hmmconcat = os.path.join(hmmdir, "all.hmm")
        # concatenate
        if len(profiles) == 0:
            self.stop("Placement resulted in zero profiles")
            return(False)
        
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
        hmmDir = os.path.join(self.cfg['outdir'], "workfiles", "hmmer", "estimations")
        file.isdir(hmmDir)
        hmmOut = os.path.join(hmmDir, "placement.tsv")
        hmmOus = os.path.join(hmmDir, "placement.out")
        hitOut = os.path.join(hmmDir, "hits.tsv")
        
        h = hmmer("hmmsearch", proteinfaa, hmmOut)
        if h.doIneedTorun(self.cfg['force']) or self.cfg['fplace'] or file.isnewer(hmmfile, hmmOut):
            log("Running hmmer for chosen locations", self.cfg['verbose'])
            h.run(hmmOus, hmmfiles = hmmfile, 
                  evalue = self.cfg['evalue'],
                  cores = self.cfg['threads'])
            # clean hmmer outpout
            log("Processing Hmmer results", self.cfg['verbose'])
            hitOut = h.clean(hmmOut, bedfile, hitOut, self.cfg['mindist'])
        return(hitOut)
         
    def gmes(self, fasta):
        """
        predict proteins using gmes
        """
        gmesDir = os.path.join(self.cfg['outdir'], "workfiles", "gmes")
        file.isdir(gmesDir)
        gmesOut = os.path.join(gmesDir, "prot_seq.faa")
        gtffile = os.path.join(gmesDir, "genemark.gtf")
        inputfasta = os.path.abspath(os.path.join(gmesDir, "input.fna"))
            
        # GeneMark-ES
        g = gmes("runGMES", fasta, gmesOut)
        if g.doIneedTorun(self.cfg['force']):
            
            if self.cfg['cleanfasta']:
                # rename fasta entries, so we dont have white spaces in them
                # can be turned of via cleanfasta in config file
                log("Copying fasta and cleaning names. Disable via 'cleanfasta' setting", self.cfg['verbose'])
                g.input = base.clearFastaNames(fasta, inputfasta)
                
            log("Running GeneMark-ES", self.cfg['verbose'])
            g.run(cores = self.cfg['threads'])
            
        # always check if gtffile exists, if not Genemark-ES failed and 
        # we can stop here
        if not file.exists(gtffile):
            log("GeneMark-ES failed on this bin")           
            self.stop("GeneMark-ES failed")
            return(False, False)        
        
        # make a bed file from GTF
        bedf = os.path.join(gmesDir, "proteins.bed")
        if self.cfg['force'] or file.isnewer(gtffile, bedf):   
            log("Extracting protein locations", self.cfg['verbose'])
            bedf = base.gmesBED(gtffile, bedf)

        return(gmesOut, bedf)      
    
    def place(self, fasta, bedfile):
        """
        main function to place a bin in the tree.
        will subsequently run hmmer 
        """
        # test if we can open the input files first
        if not base.exists(fasta):
            self.stop("Could not open fasta file")
            return(False)
        if not base.exists(bedfile):
            self.stop("Could not open bedfile file")
            return(False)
        
        
        # define output files
        hmmDir = os.path.join(self.cfg['outdir'], "workfiles", "hmmer")
        file.isdir(hmmDir)
        hmmOut = os.path.join(hmmDir, "placement.tsv")
        hmmOus = os.path.join(hmmDir, "placement.out")
        hitOut = os.path.join(hmmDir, "hits.tsv")
        
        # run hmmer if forced or input newer than output            
        h = hmmer("hmmsearch", fasta, hmmOut)
        if h.doIneedTorun(self.cfg['force']) or self.cfg['fplace']:
            log("Searching for proteins to place in the tree (HMMER)", self.cfg['verbose'])
            h.run(hmmOus, hmmfiles = self.config.placementHMMs,
                  evalue = self.cfg['evalue'],
                  cores = self.cfg['threads'])
            # clean hmmer outpout
            log("Processing Hmmer results", self.cfg['verbose'])
            hitOut = h.clean(hmmOut, bedfile, hitOut, self.cfg['mindist'])
       
        
        # pplacer paths
        placerDir = os.path.join(self.cfg['outdir'],"workfiles","pplacer")
        placerDirTmp = os.path.join(placerDir, "tmp")
        pplaceAlinment = os.path.join(placerDir, "horizontalAlignment.fasta")
        pplaceOut = os.path.join(placerDir, "placement.jplace")
        pplaceOutReduced = os.path.join(placerDir, "placementReduced.jplace")
        file.isdir(placerDirTmp)
        
        # pplacer
        log("Preparing pplacer", self.cfg['verbose'])
        pp = pplacer("pplacer", fasta, pplaceOut)
        if pp.doIneedTorun(self.cfg['force']) or self.cfg['fplace']:
            log("Creating alignments", self.cfg['verbose'])
            pp.prepareAlignment(pplaceAlinment, hitOut, os.path.join(self.config.dirname, "profile.list"), fasta, 
                                self.config, self.cfg,  placerDirTmp )
            if pp.lenscmgs == 0:
                self.stop("No protein sequences to place")
            else:
                log("Placing alignments", self.cfg['verbose'])
                pp.run(os.path.join(self.config.dirname, "refpkg", "concat.refpkg"),
                       cores = self.cfg['threads'])
        
        # reduce placements to the placements with at least posterior of p
        log("reducing placements using likelyhood", self.cfg['verbose'])
        pplaceOutReduced = pp.reduceJplace(pplaceOut, pplaceOutReduced, self.cfg['minPlacementLikelyhood'])
        
        # run TOG to get a tree
        togTree = os.path.join(placerDir, "placement.tree")
        tg = tog("guppy", pplaceOutReduced, togTree)
        if tg.doIneedTorun(self.cfg['force']):
            log("Fetching pplacer tree (guppy tog)", self.cfg['verbose'])
            tg.run()
        
        log("Getting the best placement(s)", self.cfg['verbose'])
        # now we can place the bin using the tree
        t = treelineage.treeHandler(togTree)
        t2 = treelineage.treeHandler(self.config.tree)
        orignialleaves = t2.leaves()
        sets = self.getSets()
        placements = t.getPlacement(self.cfg['placementMethod'], 
                                    sets, 
                                    orignialleaves, 
                                    self.cfg['nPlacements'], 
                                    self.cfg['minSupport'])
        log("Done placing, continuing with quality estimates", self.cfg['verbose'])
        return(placements)
        
    def getSets(self):
        setinfo = os.path.join(self.config.dirname, "sets", "setinfo.csv")
        # load sets and reduce to sets matching parameters 
        sets = []
        ints = ["n", "ngenomes"]
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
                    # convert int columns to int instead of having them as str
                    for k in ints:
                        n[k] = int(n[k])
                    sets.append(n)
                    
        return(sets)
                
        
        
        
