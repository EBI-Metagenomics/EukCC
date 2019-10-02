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
import hashlib

from ete3 import NCBITaxa
ncbi = NCBITaxa()

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
                 debug = None,
                 place = None, 
                 verbose = True,
                 evalue = None,
                 force = None, 
                 fplace = None,
                 cleanfasta = None,
                 isprotein = None, 
                 bedfile = None,
                 hmm = None,
                 noplace = None,
                 training = None
                ):
        # check config dir
        self.config = eukinfo(configdir)
        self.cfg = self.config.cfg
        
        # update config with function params
        self.cfg = updateConf(self.cfg, "cleanfasta", cleanfasta)
        self.cfg = updateConf(self.cfg, "evalue", evalue)
        self.cfg = updateConf(self.cfg, "force", force)
        self.cfg = updateConf(self.cfg, "fplace", fplace)
        self.cfg = updateConf(self.cfg, "isprotein", isprotein)
        self.cfg = updateConf(self.cfg, "hmm", hmm)
        self.cfg = updateConf(self.cfg, "noplace", noplace)
        self.cfg = updateConf(self.cfg, "debug", debug)
        self.cfg = updateConf(self.cfg, "outdir", outdir)
        self.cfg = updateConf(self.cfg, "outfile", outfile)
        self.cfg = updateConf(self.cfg, "place", place)
        self.cfg = updateConf(self.cfg, "training", training)
        self.cfg = updateConf(self.cfg, "threads", threads)
        self.cfg = updateConf(self.cfg, "verbose", verbose)
        
        self.stopped = {"stopped": False,
                        "reason": ""}
        
        # set name of run
        self.cfg['name'] = name = (os.path.splitext(os.path.basename(fastapath))[0])
        
        # if training, we need so set E-Value to trainingEvalue
        if self.cfg['training']:
            self.cfg['evalue'] = self.cfg['trainingEvalue']

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
        if (self.cfg['training'] or self.cfg['hmm']) and not self.stopnow():
            log("Running on custom hmm", self.cfg['verbose'])
            _a = self.runPlacedHMM(self.cfg['hmm'], proteinfaa, bedfile)
            
        
        if self.cfg['noplace']:
            self.stop("Stopping because we dont want to place the sequence")
        
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
            # infer lineage
            _ = self.inferLineage(self.placed[self.cfg['placementMethod']])
            _ = self.plot()
            
        # estimate completeness and contamiantion
        if not self.stopnow():            
            outputfile = os.path.join(self.cfg['outdir'], self.cfg['outfile'])
            self.estimate(hits, outputfile, self.placed[self.cfg['placementMethod']])
    
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
    
    def updateStep(self, step, status):
        self.config.cfg[step] = status
    
    def writeProcess(self):
        print(self.config.cfg['GeneMark-ES'])
    
    def estimate(self, hits, outfile, placements):
        hit = {}
        log("Estimating scores now")
        
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
            s = self.readSet(placements[i]['node'])
            # completeness is the overap of both sets 
            cmpl = len(singletons & s)/len(s)
            cont = len(multitons & s)/len(s)

            # make to percentage and round to 2 positions
            placements[i]['completeness'] = round(cmpl * 100, 2)
            placements[i]['contamination'] = round(cont * 100, 2)
        
        log("finished estimating", self.cfg['verbose'])
        
        # write to output file
        k = ["completeness", "contamination", "node", "n", 
             "ngenomes", "cover", "nPlacements",
             "taxid", "lineage", "taxidlineage", "file"]
        with open(outfile, "w") as f:
            f.write("{}\n".format("\t".join(k)))
            for p in placements:
                # insert the file name
                p['file'] = self.cfg['name']
                # write to file
                f.write("{}\n".format("\t".join([str(p[key]) for key in k])))
                
        log("Wrote estimates to: {}".format(outfile), self.cfg['verbose'])
        
        # done
        return(True)
    
    def readSet(self, node):
        profiles = []
        localpath = os.path.join(self.config.dirname, "sets", "{}.set".format(node))
        with open(localpath) as f:
            for line in f:
                profiles.append(line.strip())
        return(set(profiles))        
    
    def concatHMM(self, places):
        profiles = []
        for p in places[self.cfg['placementMethod']]:
            localpath = os.path.join(self.config.dirname, "sets", "{}.set".format(p['node']))
            with open(localpath) as f:
                for line in f:
                    profiles.append(line.strip())
       
        
        # create all paths for all hmms
        hmmerpaths = [os.path.join(self.config.dirname, "hmms", "panther", "{}.hmm".format(profile)) 
                      for profile in profiles]
        
        # create a dir for this
        hmmdir = os.path.join(self.cfg['outdir'],"workfiles","hmmer", "estimations")
        file.isdir(hmmdir)
        hmmconcat = os.path.join(hmmdir, "all.hmm")
        
         # sort and check if we already have the hmm for this
        profiles.sort()
        canuseprev = False
        profilehash = hashlib.sha256("_".join(profiles).encode()).hexdigest()
        hashpath = os.path.join(hmmdir, "all.hash")
        if file.exists(hashpath):
            with open(hashpath) as f:
                for line in f:
                    prevhash = line.strip()
                    break
            canuseprev = prevhash == profilehash
        
        if canuseprev:
            # we can use the existing file, so no need to continue
            log("Using pressed hmms from last run")
            return(hmmconcat)
        
        # concatenate
        if len(profiles) == 0:
            self.stop("Placement resulted in zero profiles")
            return(False)
        
        log("{} hmm profiles need to be used for estimations".format(len(profiles)), self.cfg['verbose'])
        log("Concatenating hmms, this might take a while (IO limited)", self.cfg['verbose'])
        hmmconcat = base.concatenate(hmmconcat, hmmerpaths)
        # press
        log("Pressing hmms", self.cfg['verbose'])
        hp = hmmpress("hmmpress", hmmconcat, None)
        hp.run()
        
        # save profile hash
        with open(hashpath, "w") as f:
            f.write(f"{profilehash}")
        
        return(hmmconcat)

    def runPlacedHMM(self, hmmfile, proteinfaa, bedfile):
        # run hmmer and strip down

        # define output files
        hmmDir = os.path.join(self.cfg['outdir'], "workfiles", "hmmer", "estimations")
        file.isdir(hmmDir)
        hmmOut = os.path.join(hmmDir, "placement.tsv")
        hmmOus = os.path.join(hmmDir, "placement.out")
        hitOut = os.path.join(hmmDir, "hits.tsv")
            
        h = hmmer("hmmsearch", proteinfaa, hmmOut, self.cfg['debug'])
        if h.doIneedTorun(self.cfg['force']) or self.cfg['fplace'] or file.isnewer(hmmfile, hmmOut):
            log("Running hmmer for chosen locations", self.cfg['verbose'])
            h.run(hmmOus, hmmfiles = hmmfile, 
                  evalue = self.cfg['evalue'],
                  cores = self.cfg['threads'],
                  training = self.cfg['training'])
            # clean hmmer outpout
            log("Processing Hmmer results", self.cfg['verbose'])
            hitOut = h.clean(hmmOut, bedfile, hitOut, self.cfg['mindist'])
        return(hitOut)
    
    def inferLineage(self, places):
        '''
        infer the lineage from looking at the location of placement
        looking at the leaves and their tax id
        and looking at the lineages of all these
        '''
       
        # fetch file and load taxinformation
        seqinfo = self.config.pkgfile("concat.refpkg", "seq_info")
        taxids = {}
        si = base.readCSV(seqinfo)
        # make dictionary 
        for r in si:
            taxids[r['seqname']] = r["tax_id"]
        # load tree
        tree = treelineage.treeHandler(self.config.tree, annotate = False)
        # for each placement:
        for p  in places:
            # get the GCA names
            children = p['sisters'] #tree.children(p['tax_id'])
            if self.cfg['debug']:
                print(f"Infered children {children}")
            # fetch lineages for all
            lngs = []
            for c in children:
                lngs.append(ncbi.get_lineage(taxids[c]))
            
            # find common elements:
            common = set(lngs[0])
            for l in lngs[1:]:
                common = common & set(l)
 
            # common lineage
            lng = []
            for v in lngs[0]:
                if v not in common:
                    break
                # add common elements
                lng.append(v)
          
            
            nodetaxid = lng[-1]
            # now we can make it pretty
            if self.cfg['lineage'] == "limited":
                # limit to desired ranks2
                desired_ranks  = ['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
                lineage2ranks = ncbi.get_rank(lng)
                ranks2lineage = dict((rank, taxid) for (taxid, rank) in lineage2ranks.items())
                ranks  = {'{}_id'.format(rank): ranks2lineage.get(rank, 'NA') for rank in desired_ranks}
                lng = [i for i in lng if i in ranks.values()]
            # get translator and make string
            named = ncbi.translate_to_names(lng)
            # save to placed object
            p['lineage'] = "_".join(named)
            # replace names with spaces into points
            p['lineage'] = p['lineage'].replace(" ", ".")
            p['taxidlineage'] = "_".join([str(x) for x in lng])
            p['taxid'] = nodetaxid
            
            if self.cfg['debug']:
                print("infered lineage")
                print(p['lineage'])
     
        return()
    
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
            self.updateStep('gmes', 'started')
            g.run(cores = self.cfg['threads'])
            
        # always check if gtffile exists, if not Genemark-ES failed and 
        # we can stop here
        if not file.exists(gtffile):
            # log and document failing
            # then stop pipeline
            log("GeneMark-ES failed on this bin")           
            self.stop("GeneMark-ES failed")
            self.updateStep('gmes', 'failed')
            return(False, False)        
        
        # make a bed file from GTF
        bedf = os.path.join(gmesDir, "proteins.bed")
        if self.cfg['force'] or file.isnewer(gtffile, bedf):   
            log("Extracting protein locations", self.cfg['verbose'])
            bedf = base.gmesBED(gtffile, bedf)
            self.updateStep('gmes', 'created bed')
        # update final step
        self.updateStep('gmes', 'finished')

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
            self.updateStep('findprots', 'looked for proteins')
        
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
            pp.prepareAlignment(pplaceAlinment, hitOut, 
                                os.path.join(self.config.dirname, "profile.list"), 
                                fasta, 
                                self.config, self.cfg,  placerDirTmp )
            if pp.lenscmgs == 0:
                self.stop("No protein sequences to place")
                self.updateStep('findprots', 'no single copy markers')
            else:
                log("Placing alignments", self.cfg['verbose'])
                self.updateStep('pplacer', 'starting')
                pp.run(os.path.join(self.config.dirname, "refpkg", "concat.refpkg"),
                       cores = self.cfg['threads'])
                
                
        # reduce placements to the placements with at least posterior of p
        log("reducing placements using likelyhood", self.cfg['verbose'])
        if self.cfg['debug']:
            print(pplaceOutReduced)
        pplaceOutReduced = pp.reduceJplace(pplaceOut, pplaceOutReduced, self.cfg['minPlacementLikelyhood'])
        if self.cfg['debug']:
            print("Done plceing reduced")
        # run TOG to get a tree
        togTree = os.path.join(placerDir, "placement.tree")
        tg = tog("guppy", pplaceOutReduced, togTree)
        if tg.doIneedTorun(self.cfg['force']):
            log("Fetching pplacer tree (guppy tog)", self.cfg['verbose'])
            tg.run()
        
        log("Getting the best placement(s)", self.cfg['verbose'])
        # save path to togtree for plotting later
        self.cfg['togtreepath'] = togTree
        self.cfg['togjson'] = pplaceOutReduced
        # now we can place the bin using the tree
        t = treelineage.treeHandler(togTree, annotate = False)
        t2 = treelineage.treeHandler(self.config.tree, annotate = False)
        #orignialleaves = t2.leaves()
        sets = self.getSets()
        # get HCA and LCA placements
        placements = {}
        for method in ['LCA', 'HPA']:
            placements[method] = t.getPlacement(method, 
                                    sets, 
                                    t2,
                                    self.cfg['nPlacements'], 
                                    self.cfg['minSupport'],
                                    self.cfg['debug'])

        self.updateStep('pplacer', 'finished')
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


    def plot(self):
        # load tree
        log("Plotting a trees of placements")
        t = treelineage.treeHandler(self.cfg['togtreepath'], annotate = False)
        t.plot(self.placed, 
               self.cfg['togjson'],
               self.cfg['outdir'],
               self.cfg)
        return()
        
