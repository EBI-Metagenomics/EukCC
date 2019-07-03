import os
import subprocess
from eukcc.fileoperations import file
from eukcc import base
import operator 
import json

from pyfaidx import Fasta

#from collections import OrderedDict

class run():
    '''
    class to handle and run external software
    check if a software exists
    handle input and output
    '''
    def __init__(self, program, inf, outf):
        # check software is in path:
        if run.which(program) is None:
            print("{} is not installed".format(program))
        self.program = program
        self.input = inf
        self.output = outf
        if outf is not None:
            # create output dir
            dc = file.isdir(os.path.dirname(self.output))

    def doIneedTorun(self, force=False):
        if force:
            return(True)
        else:
            return(file.isnewer(self.input, self.output))

    def which(program):
        # taken from
        # https://stackoverflow.com/questions/377017/\
        # test-if-executable-exists-in-python
        def is_exe(fpath):
            return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

        fpath, fname = os.path.split(program)
        if fpath:
            if is_exe(program):
                return program
        else:
            for path in os.environ["PATH"].split(os.pathsep):
                exe_file = os.path.join(path, program)
                if is_exe(exe_file):
                    return exe_file

        return None

    def run(self, inp, outp, params=[]):
        lst = [self.program]
        lst.extend(params)
        lst.append(inp)
        try:
            with open(outp, "w") as f:
                subprocess.run(lst, stdout=f, check=True)
            return(True)
        except subprocess.CalledProcessError as e:
            print("an error occured while executing {}".format(self.program))
            print(e.output)
            return(False)
        
    def checkOutput(self, file, comment = "#"):
        """
        generic function to see if we have at least one none comment line
        in the output
        """
        with open(file) as f:
            for line in f:
                if not line.startswith(comment):
                    return(True)
        return(False)


# defining classes based on run for executing gmes and hmmscan
class gmes(run):
    def run(self, cores=1):
        lst = [self.program, "-i", self.input,
               "-o", self.output,
               "-n", str(cores)]
        try:
            subprocess.run(lst,  check=True, shell=False)
            return(True)
        except subprocess.CalledProcessError:
            print("an error occured while executing {}".format(self.program))
            return(False)

        
class hmmpress(run):
    def run(self):
        lst = [self.program, "-f", self.input]
        try:
            subprocess.run(lst,  check=True, shell=False)
            return(True)
        except subprocess.CalledProcessError:
            print("an error occured while executing {}".format(self.program))
            return(False)

class hmmer(run):
    def run(self, stdoutfile, hmmfiles, cores=1, evalue = 1e-5):
        lst = [self.program, "--cpu", str(cores),
               "-E", str(evalue),
               "-o", stdoutfile,
               "--tblout", self.output,
               hmmfiles, self.input]
        #print(" ".join(lst))
        try:
            subprocess.run(lst,  check=True, shell=False)
            return(True)
        except subprocess.CalledProcessError:
            print("an error occured while executing {}".format(self.program))
            return(False)
    

    def clean(self, hmmer, bedfile, resultfile, mindist = 2000):
        bed = base.readbed(bedfile)
        fields = {"subject": 0,
                  "evalue": 4,
                  "profile": 2}
        table = []
        # intcoluns = 
        ints  = ['start', 'stop']
        # read in Hmmer results
        with open(hmmer) as h:
            for line in h:
                if line.startswith("#"):
                    continue
                l = line.split()
                n = {}
                for k,v in fields.items():
                    n[k] = l[v]
                for k in ['chrom', 'start', 'stop']:
                    if n['subject'] in bed.keys():
                        n[k] = bed[n['subject']][k]
                    else:
                        print("Could not find entry in bed file for {}".format(n['subject']))
                # fix profile name
                n['profile'] = n['profile'].split(".")[0]
                n['keep'] = True
                
                # convert int columns to int instead of having them as str
                for k in ints:
                    n[k] = int(n[k])
                
                table.append(n)
        
        # sort table by profile, chromsome, start, end
        columns = ["profile", "chrom", "start", "stop"]
        columns.reverse()
        for k in columns:
            table.sort(key=operator.itemgetter(k))
        
        def closeOrWithin(row, lastrow, distanceTresh):
            '''
            returns true if row is within or close to the lastrow
            '''
            # start within
            if row["start"] >= lastrow["start"] and row["start"] <= lastrow["stop"]:
                return(True)
            # end within
            if row["stop"] >= lastrow["start"] and row["stop"] <= lastrow["stop"]:
                return(True)

            # distance between any value with any other value:
            a = [row["start"], row["stop"]]
            b = [lastrow["start"], lastrow["stop"]]
            for v in a:
                for w in b:
                    if abs(v - w) <= distanceTresh:
                        return(True)
    
            return(False)
        
        # now we can iterate this list and alway retain the one with the higher evalue        
        j = 0
        for i in range(1, len(table)):
            # set values for this iteration
            row = table[i]
            lastrow = table[j]
            keepI = True
            keepJ = True
            # check if hits are of same profile and on same chromsome
            # only then we can merge
            if row['profile'] == lastrow['profile'] and row['chrom'] == lastrow['chrom']:
                # check if both are close
                if closeOrWithin(row, lastrow, mindist):
                    if float(lastrow['evalue']) < float(row['evalue']):
                        keepI = False
                    else:
                        keepJ = False
                
            
            # save which row we keep
            if "keep" not in table[j].keys() or keepJ == False:
                table[j]['keep'] = keepJ
            if "keep" not in table[i].keys() or keepI == False:
                table[i]['keep'] = keepI
            
            # update j, so index to lastrow
            if keepJ == False:
                j = i 
            elif keepJ == True and keepI == True:
                j = i
            elif keepJ == True and keepI == False:
                # j stays the same
                j = j # obv redundant, but nice for making sense of this code
        
            
        # write result to file
        cols = ['profile', 'subject', 'chrom', 'start', 'stop', 'evalue']
        with open(resultfile, "w") as f:
            f.write("\t".join(cols + ["\n"]))
            for row in table:
                if row['keep'] == False:
                    continue
                # print if we keep the row
                v = [str(row[k]) for k in cols] + ["\n"]
                f.write("\t".join(v))
        
        return(resultfile)
    

    

class hmmalign(run):
    def run(self, alignmentpath, hmmpath):
        lst = [self.program, "--outformat", "afa", 
               "--mapali",
               alignmentpath,  
               "--amino", hmmpath, self.input]
        #print(" ".join(lst))
        try:
            with open(self.output, "w") as f:
                subprocess.run(lst,  check=True, shell=False, stdout = f)
            return(True)
        except subprocess.CalledProcessError:
            print("an error occured while executing {}".format(self.program))
            return(False)

    
class pplacer(run):
    def run(self, pkg, cores = 1):
        lst = [self.program, 
                    "-o", self.output ,"-p", "--keep-at-most", "5",
                    "-m", "LG",
                    "-j", str(cores), 
                    "-c", pkg , self.input]
        print(" ".join([str(i) for i in lst]))
        try:
            subprocess.run(lst,  check=True, shell=False)
            return(True)
        except subprocess.CalledProcessError:
            print("an error occured while executing {}".format(self.program))
            return(False)
        
    def prepareAlignment(self, pplaceAlinment, hmmerOutput, proteinList, proteinFasta, config, cfg, tmpDir ):
        # for each profile that we found a SCMG
        # make concatenated fasta (including new sequence)
        # align that to the hmm profile
        profiles = []
        with open(proteinList) as pl:
            for l in pl:
                profiles.append(l.strip())
        # read in hmmer output to get seqnames of target proteins
        cols = []
        #hmmer = []
        #scmgsr = []

        #with open(hmmerOutput) as ho:
        #    for line in ho:
        #        l = line.split()
        #        if len(cols) == 0:
        #            cols = l
        #        else:
        #            n = {}
        #            for k,v in zip(cols, l):
        #                n[k] = v
        #            hmmer.append(n)
        #            scmgsr.append(n['profile'])
        
        # load hmmer results
        hmmer = base.readTSV(hmmerOutput)
        # extract SCMGs
        scmgsr = [hit['profile'] for hit in hmmer]
        scmgs = [p for p in set(scmgsr) if scmgsr.count(p) == 1]
        scmgs.sort()
        # count the placements
        self.lenscmgs = len(scmgs)
        
        if self.lenscmgs > 0:
            # get the protein names for each scmgs
            proteinnames = []
            for profile in scmgs:
                for row in hmmer:
                    if row['profile'] == profile:
                        proteinnames.append(row['subject'])
            # load proteins
            proteins = Fasta(proteinFasta)

            alignments = []
            # for each protein/SCMG we need to make an alignment
            for p, g in zip(scmgs, proteinnames):
                # write seq to file
                seq = proteins[g]
                name = "{}_{}".format(p, g)
                tmpfasta = os.path.join(tmpDir, "{}.faa".format(name))
                geneAlignment = os.path.join(tmpDir, "{}.aln".format(name))
                tmpfasta = base.writeFasta(tmpfasta, name, seq)
                profileAlignment = config.pkgfile("{}.refpkg".format(p), "aln_fasta")
                profileHMM = config.pkgfile("{}.refpkg".format(p), "profile")
                # start and run alignment of this profile
                ha = hmmalign("hmmalign", tmpfasta, geneAlignment)
                ha.run(profileAlignment, profileHMM)
                alignments.append(geneAlignment)        

            # after this we concatenate the alignments
            # ordered by the profile name
            self.input = base.horizontalConcat(pplaceAlinment, alignments, scmgs, 
                             config.pkgfile("{}.refpkg".format("concat"), "aln_fasta"))
    
    def reduceJplace(self, jplace, jplaceselection, placementCutoff = 0.5):
        '''
        reduces a jplace file to placements with at least placementCutoff
        post_prob
        '''        
        with open(jplace) as json_file:  
            j = (json.load(json_file))
            fields = j['fields']
            newplacements = []
            for placement in j['placements']:
                ps = placement['p']
                kp = []
                for p in ps:
                    d = {k: v for k,v in zip(fields, p)}
                    if d['post_prob'] >= placementCutoff:
                        kp.append(p)

                if len(kp) > 0:
                    placementn = placement.copy()
                    placementn['p'] = kp
                    newplacements.append(placementn)
            jn = j.copy()
            jn['placements'] = newplacements

        with open(jplaceselection, "w") as f:
            json.dump(jn,f)
        return(jplaceselection)
        

        
class tog(run):
    def run(self):
        lst = [self.program, "tog", self.input, "-o", self.output]
        try:
            subprocess.run(lst,  check=True, shell=False)
            return(True)
        except subprocess.CalledProcessError:
            print("an error occured while executing {}".format(self.program))
            return(False)

        
def checkDependencies(dep):
    """
    simple function to check all dependencies
    """
    allgood = True
    for name, program in dep.items():
        if run.which(program) is None:
            print("Dependency: Could not find {} ({})".format(program, name))
            allgood = False

    return(allgood)
