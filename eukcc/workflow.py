#!/usr/bin/env python3
#
# This file is part of the EukCC (https://github.com/openpaul/eukcc).
# Copyright (c) 2019 Paul Saary
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
# provides all file operation functions
# used inthis package
import os
import logging
from eukcc.info import eukinfo
from eukcc.base import log
from eukcc.exec import hmmer
from eukcc.exec import hmmpress
from eukcc.exec import gmes
from eukcc.exec import pplacer
from eukcc.exec import tog
from eukcc import base
from eukcc import treelineage
from eukcc.fileoperations import file
import hashlib
from ete3 import NCBITaxa
from pyfaidx import Fasta
from collections import defaultdict


dep = {
    "hmmer": "hmmsearch",
    "runGMES": "runGMES",
    "GeneMarkES": "gmes_petap.pl",
    "pplacer": "pplacer",
}


class eukcc:
    def __init__(self, options):
        # define config class
        # stores all magic numbers and user settings
        self.config = eukinfo(options)
        self.cfg = self.config.cfg
        self._clean_fasta = None

        # set name of run
        self.cfg["name"] = os.path.splitext(os.path.basename(options.fasta))[0]
        if self.cfg["training"]:
            self.cfg["evalue"] = self.cfg["trainingEvalue"]
            logging.info("Defining e-value for training as %", self.cfg["evalue"])
        elif self.cfg["dbinfo"]["modus"] == "evalue":
            self.cfg["evalue"] = float(self.cfg["dbinfo"]["evalue"])
            logging.debug("Defining e-value to match DB info: %f", self.cfg["evalue"])

    def checkIO(self, fastapath, outdir):
        # create outdir if not exists
        file.isdir(self.cfg["outdir"])
        # check if input and output can be accessed
        logging.debug("Warning: IO check not yet implemented")
        return False

    def updateStep(self, step, status):
        self.config.cfg[step] = status

    def writeProcess(self):
        print(self.config.cfg["GeneMark-ES"])

    def get_silent_contig(self, fasta, hits, placements):
        """
        given hits and the fasta file, we can calculate how many
        DNA bp can not be detected as contaminant, because no markers
        are on these contigs
        """
        logging.debug("Calculating silent contig fraction")
        if fasta is None:
            logging.debug("As no DNA fasta was given, we can't compute the contig fraction")
            for placement in placements:
                placement["max_silent_contamination"] = "NA"
            return

        # get dict to link profile to contigs
        profile_2_contigs = defaultdict(set)
        for row in base.readTSV(hits):
            profile_2_contigs[row["profile"]].add(row["chrom"])

        # read in fasta once
        dna_lens = {}
        for rec in Fasta(fasta):
            dna_lens[rec.name] = len(str(rec))
        total_len = sum([v for k, v in dna_lens.items()])

        # for each placement compute the silent fraction
        for placement in placements:
            contigs = set()
            # make union of contigs
            for profile in placement["set"]:
                contigs = contigs | profile_2_contigs[profile]
            # define missing contigs as contig names we did not locate any marker genes on
            m_contigs = set(dna_lens.keys()) - contigs
            missing_len = sum([dna_lens[contig] for contig in m_contigs])
            placement["max_silent_contamination"] = round(missing_len / total_len * 100, 2)
            logging.debug("The silent fraction could be up to {}%".format(placement["max_silent_contamination"]))
        return

    def estimate(self, hits, outfile, placements):
        hit = {}
        logging.info("Estimating scores now")

        if self.cfg["touch"]:
            file.touch(outfile)
            logging.info("Returning as we only touch")
            return

        r = base.readTSV(hits)
        # count profile hits
        for row in r:
            if row["profile"] not in hit.keys():
                hit[row["profile"]] = 1
            else:
                hit[row["profile"]] += 1

        singletons = set(hit.keys())
        multitons = set([k for k, v in hit.items() if v > 1])

        # now we can estimate completeness and contamination for each placement
        for i in range(len(placements)):
            s = self.readSet(placements[i]["node"])
            placements[i]["set"] = s
            # completeness is the overap of both sets
            cmpl = len(singletons & s) / len(s)
            cont = len(multitons & s) / len(s)

            # make to percentage and round to 2 positions
            placements[i]["completeness"] = round(cmpl * 100, 2)
            placements[i]["contamination"] = round(cont * 100, 2)

        # compute silent fraction per placement and set
        self.get_silent_contig(self._clean_fasta, hits, placements)

        log("Finished estimating")
        self.write_outfile(outfile, placements)

        # done
        return True

    def write_outfile(self, outfile=None, result=None):
        if outfile is None:
            outfile = os.path.join(self.cfg["outdir"], "eukcc.tsv")
        # write to output file
        k = [
            "completeness",
            "contamination",
            "max_silent_contamination",
            "node",
            "n",
            "ngenomes",
            "cover",
            "nPlacements",
            "taxid",
            "lineage",
            "taxidlineage",
            "file",
        ]
        with open(outfile, "w") as f:
            f.write("{}\n".format("\t".join(k)))
            if result is None:
                logging.warning("No estimates were written")
                exit(11)
            for p in result:
                # insert the file name
                p["file"] = self.cfg["name"]
                # write to file
                f.write("{}\n".format("\t".join([str(p[key]) for key in k])))

        log("Wrote estimates to: {}".format(outfile))

    def readSet(self, node):
        profiles = []
        localpath = os.path.join(self.cfg["db"], "sets", "{}.set".format(node))
        with open(localpath) as f:
            for line in f:
                profiles.append(line.strip())
        return set(profiles)

    def concatHMM(self):
        # create a dir for this
        hmmdir = os.path.join(self.cfg["outdir"], "workfiles", "hmmer", "estimations")
        file.isdir(hmmdir)
        hmmconcat = os.path.join(hmmdir, "all.hmm")

        if self.cfg["touch"]:
            file.touch(hmmconcat)
            return hmmconcat

        profiles = set()
        for p in self.placements[self.cfg["placementMethod"]]:
            localpath = os.path.join(self.cfg["db"], "sets", "{}.set".format(p["node"]))
            with open(localpath) as f:
                for line in f:
                    profiles.add(line.strip())
        # make profiles to sorted list
        profiles = list(profiles)
        profiles.sort()
        # create all paths for all hmms
        hmmerpaths = [os.path.join(self.cfg["db"], "hmms", "panther", "{}.hmm".format(profile)) for profile in profiles]
        # sort and check if we already have the hmm for this
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
            logging.info("Using pressed hmms from last run")
            return hmmconcat

        # concatenate
        if len(profiles) == 0:
            logging.error("We have no profiles to evaluate")
            self.write_outfile()
            exit(1)

        log("{} hmm profiles need to be used for estimations".format(len(profiles)))
        log("Concatenating hmms, this might take a while (IO limited)")
        hmmconcat = base.concatenate(hmmconcat, hmmerpaths)
        # press
        log("Pressing hmms")
        hp = hmmpress("hmmpress", hmmconcat, None, touch=self.cfg["touch"])
        hp.run()

        # save profile hash
        with open(hashpath, "w") as f:
            f.write(f"{profilehash}")

        return hmmconcat

    def runPlacedHMM(self, hmmfile, proteinfaa, bedfile):
        # run hmmer and strip down

        # define output files
        hmmDir = os.path.join(self.cfg["outdir"], "workfiles", "hmmer", "estimations")
        file.isdir(hmmDir)
        hmmOut = os.path.join(hmmDir, "placement.tsv")
        hmmOus = os.path.join(hmmDir, "placement.out")
        hitOut = os.path.join(hmmDir, "hits.tsv")

        h = hmmer("hmmsearch", proteinfaa, hmmOut, self.cfg["debug"], touch=self.cfg["touch"],)
        if h.doIneedTorun(self.cfg["force"]) or self.cfg["fplace"] or file.isnewer(hmmfile, hmmOut):
            logging.info("Running hmmer for chosen locations")
            h.run(
                hmmOus,
                hmmfiles=hmmfile,
                modus=self.cfg["dbinfo"]["modus"],
                evalue=self.cfg["evalue"],
                cores=self.cfg["ncores"],
                training=self.cfg["training"],
            )
            # clean hmmer outpout
            logging.info("Processing Hmmer results")
            hitOut = h.clean(hmmOut, bedfile, hitOut, self.cfg["mindist"])
        return hitOut

    def inferLineage(self, places):
        """
        infer the lineage from looking at the location of placement
        looking at the leaves and their tax id
        and looking at the lineages of all these
        """
        if self.cfg["touch"]:
            return
        ncbi = NCBITaxa()
        # fetch file and load taxinformation
        seqinfo = self.config.pkgfile("concat.refpkg", "seq_info")
        taxids = {}
        si = base.readCSV(seqinfo)
        # make dictionary
        for r in si:
            taxids[r["seqname"]] = r["tax_id"]

        # for each placement:
        logging.debug("Infering lineages now")
        for p in places:
            # get the GCA names
            children = p["sisters"]
            # fetch lineages for all
            lngs = []
            for c in children:
                try:
                    lngs.append(ncbi.get_lineage(taxids[c]))
                except ValueError as e:
                    logging.warning(e)

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
            if not self.cfg["fullineage"]:
                # limit to desired ranks2
                desired_ranks = [
                    "superkingdom",
                    "kingdom",
                    "phylum",
                    "class",
                    "order",
                    "family",
                    "genus",
                    "species",
                ]
                lineage2ranks = ncbi.get_rank(lng)
                ranks2lineage = dict((rank, taxid) for (taxid, rank) in lineage2ranks.items())
                ranks = {"{}_id".format(rank): ranks2lineage.get(rank, "NA") for rank in desired_ranks}
                lng = [i for i in lng if i in ranks.values()]
            # get translator and make string
            named = ncbi.translate_to_names(lng)
            # save to placed object
            p["lineage"] = "_".join(named)
            # replace names with spaces into points
            p["lineage"] = p["lineage"].replace(" ", ".")
            p["taxidlineage"] = "_".join([str(x) for x in lng])
            p["taxid"] = nodetaxid

        return ()

    def gmes(self, fasta):
        """
        predict proteins using gmes
        """
        logging.debug("Starting gmes function")

        gmesDir = os.path.join(self.cfg["outdir"], "workfiles", "gmes")
        file.isdir(gmesDir)
        gmesOut = os.path.join(gmesDir, "prot_seq.faa")
        gtffile = os.path.join(gmesDir, "genemark.gtf")
        inputfasta = os.path.abspath(os.path.join(gmesDir, "input.fna"))

        # GeneMark-ES
        g = gmes("runGMES", fasta, [gtffile, gmesOut], touch=self.cfg["touch"])
        logging.debug("Defined gmes run")
        if g.doIneedTorun(self.cfg["force"]):
            # rename fasta entries, so we dont have white spaces in them
            # can be turned of via cleanfasta in config file
            if not self.cfg["touch"]:
                g.input = base.clearFastaNames(fasta, inputfasta)
            else:
                g.input = inputfasta

            logging.info("Running GeneMark-ES")
            g.run(cores=self.cfg["ncores"])
        else:
            logging.debug("I do not need to run gmes, output exists:")
            logging.debug(gtffile)

        # always check if gtffile exists, if not Genemark-ES failed and
        # we can stop here
        if not file.exists(gtffile):
            # log and document failing
            # then stop pipeline
            logging.error("GeneMark-ES failed on this bin")
            self.write_outfile()
            exit(1)
        elif self.cfg["clean"]:
            # clean temp dirs
            _tmpdirs = ["data", "run", "info", "output/data", "output/gmhmm"]
            tempdirs = [os.path.join(gmesDir, x) for x in _tmpdirs]
            g.cleanup(tempdirs)

        # make a bed file from GTF
        bedf = os.path.join(gmesDir, "proteins.bed")
        if self.cfg["force"] or file.isnewer(gtffile, bedf) and not self.cfg["touch"]:
            logging.info("Extracting protein locations")
            bedf = base.gmesBED(gtffile, bedf)

        # touch files expected for next step
        if self.cfg["touch"]:
            g.touch([bedf, gmesOut])
        self._clean_fasta = inputfasta
        return (gmesOut, bedf)

    def pygmes(self, fasta, db):
        outdir = os.path.join(self.cfg["outdir"], "workfiles", "pygmes")
        faafile = os.path.join(outdir, "predicted_proteins.faa")
        bedfile = os.path.join(outdir, "predicted_proteins.bed")
        self._clean_fasta = os.path.join(outdir, "gmesclean_{}".format(os.path.basename(fasta)))
        from pygmes import pygmes

        # check if we need to launch
        need_run = False
        if not file().exists(faafile) or not file().exists(bedfile):
            need_run = True
        elif file.isnewer(fasta, faafile):
            need_run = True
        if need_run:
            pygmes(fasta, outdir, db=db, clean=True, ncores=self.cfg["ncores"])
        # check if pg worked
        if os.path.exists(faafile) and os.path.exists(bedfile):
            if os.stat(faafile).st_size == 0 or os.stat(bedfile).st_size == 0:
                logging.warning("No predicted proteins")
                self.write_outfile()
                exit(1)
            else:
                return (faafile, bedfile)
        else:
            logging.warning("No predicted proteins, pyfaidx failed")
            self.write_outfile()
            exit(1)

    def place(self, fasta, bedfile):
        """
        main function to place a bin in the tree.
        will subsequently run hmmer
        """
        # test if we can open the input files first
        if not base.exists(fasta):
            logging.error("Could not open fasta file")
            self.write_outfile()
            exit(1)
        if not base.exists(bedfile):
            logging.error("Could not open bed file")
            self.write_outfile()
            exit(1)

        # define output files
        hmmDir = os.path.join(self.cfg["outdir"], "workfiles", "hmmer")
        file.isdir(hmmDir)
        hmmOut = os.path.join(hmmDir, "placement.tsv")
        hmmOus = os.path.join(hmmDir, "placement.out")
        hitOut = os.path.join(hmmDir, "hits.tsv")

        # run hmmer if forced or input newer than output
        h = hmmer("hmmsearch", fasta, hmmOut, touch=self.cfg["touch"])
        if h.doIneedTorun(self.cfg["force"]) or self.cfg["fplace"]:
            logging.info("Searching for proteins to place in the tree")
            h.run(
                hmmOus,
                hmmfiles=self.config.placementHMMs,
                modus=self.cfg["dbinfo"]["modus"],
                evalue=self.cfg["evalue"],
                cores=self.cfg["ncores"],
            )
            # clean hmmer outpout
            logging.info("Processing Hmmer results")
            hitOut = h.clean(hmmOut, bedfile, hitOut, self.cfg["mindist"])
            self.updateStep("findprots", "looked for proteins")

        # pplacer paths
        placerDir = os.path.join(self.cfg["outdir"], "workfiles", "pplacer")
        placerDirTmp = os.path.join(placerDir, "tmp")
        pplaceAlinment = os.path.join(placerDir, "horizontalAlignment.fasta")
        pplaceOut = os.path.join(placerDir, "placement.jplace")
        pplaceLog = os.path.join(placerDir, "placement.log")
        pplaceOutReduced = os.path.join(placerDir, "placementReduced.jplace")
        file.isdir(placerDirTmp)

        # pplacer
        logging.debug("Preparing pplacer")
        pp = pplacer("pplacer", fasta, pplaceOut, touch=self.cfg["touch"])
        if pp.doIneedTorun(self.cfg["force"]) or self.cfg["fplace"]:
            logging.debug("Preparing alignments")
            pp.prepareAlignment(
                pplaceAlinment,
                hitOut,
                os.path.join(self.cfg["db"], "profile.list"),
                fasta,
                self.config,
                self.cfg,
                placerDirTmp,
            )
            if pp.lenscmgs == 0 and not self.cfg["touch"]:
                logging.error("Could not find any marker genes")
                self.write_outfile()
                exit(1)
            else:
                logging.info("Placing proteins in tree")
                self.updateStep("pplacer", "starting")
                pplacer_success = pp.run(
                    os.path.join(self.cfg["db"], "refpkg", "concat.refpkg"),
                    logfile=pplaceLog,
                    cores=self.cfg["ncorespplacer"],
                )
                if pplacer_success is False:
                    logging.warning("Pplacer could not finish. Exiting now")
                    self.write_outfile()
                    exit(1)

        # reduce placements to the placements with at least posterior of p
        logging.debug("Reducing placements")
        if not self.cfg["touch"]:
            pplaceOutReduced = pp.reduceJplace(pplaceOut, pplaceOutReduced, self.cfg["minPlacementLikelyhood"])
        else:
            pp.touch([pplaceOutReduced])
        logging.debug("Reducing placements done")
        # run TOG to get a tree
        togTree = os.path.join(placerDir, "placement.tree")
        tg = tog("guppy", pplaceOutReduced, togTree, touch=self.cfg["touch"])
        if tg.doIneedTorun(self.cfg["force"]):
            logging.debug("Fetching tree")
            r = tg.run()
            if r is False:
                logging.debug("No placement found")
                self.write_outfile()

        logging.debug("Getting best placements")
        # save path to togtree for plotting later
        self.cfg["togtreepath"] = togTree
        self.cfg["togjson"] = pplaceOutReduced
        # now we can place the bin using the tree
        if not self.cfg["touch"]:
            t = treelineage.treeHandler(togTree, annotate=False)
            t2 = treelineage.treeHandler(self.config.tree, annotate=False)
            sets = self.getSets()
            # get HCA and LCA placements
            self.placements = {}
            for method in ["LCA", "HPA"]:
                self.placements[method] = t.getPlacement(
                    method,
                    sets,
                    t2,
                    self.cfg["nPlacements"],
                    self.cfg["minSupport"],
                    maximum=self.cfg["nEvals"],
                    debug=self.cfg["debug"],
                )
        else:
            self.placements = {"LCA": "touch", "HCA": "touch"}

        logging.info("MAG succesfully placed in tree")

    def getSets(self):
        setinfo = os.path.join(self.cfg["db"], "sets", "setinfo.csv")
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
                if int(n["ngenomes"]) >= self.cfg["minGenomes"] and int(n["n"]) >= self.cfg["minProfiles"]:
                    # convert int columns to int instead of having them as str
                    for k in ints:
                        n[k] = int(n[k])
                    sets.append(n)

        return sets

    def plot(self):
        # load tree
        logging.info("Plotting trees of placements")
        if self.cfg["touch"]:
            return

        t = treelineage.treeHandler(self.cfg["togtreepath"], annotate=False)
        t.plot(self.placements, self.cfg["togjson"], self.cfg["outdir"], self.cfg)
