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
import os
import subprocess
import shutil
from eukcc.fileoperations import file
from eukcc import base
import operator
import json
import logging

from pyfaidx import Fasta

# from collections import OrderedDict


class run:
    """
    class to handle and run external software
    check if a software exists
    handle input and output
    """

    def __init__(self, program, inf, outf, debug=False, touch=False):
        # check software is in path:
        if run.which(program) is None:
            print("{} is not installed".format(program))
        self.debug = debug
        self.touchonly = touch
        self.program = program
        self.input = inf
        # in case multiple output fiules are defined
        # we set the first one as output but use all for testing
        # is a rule has to be run
        if isinstance(outf, list):
            self.output = outf[0]
            self.output_test = outf
        else:
            self.output = outf
            self.output_test = [outf]

        if outf is not None:
            # create output dir
            file.isdir(os.path.dirname(self.output))

    def touch(self, files=None):
        """
        simple function, to only touch output files
        """
        logging.debug("Only touching now")
        if files is None:
            path = self.output
            with open(path, "a"):
                os.utime(path, None)
        else:
            for path in files:
                with open(path, "a"):
                    os.utime(path, None)

    def doIneedTorun(self, force=False):
        logging.debug("Testing if I need to run this step")
        if force or self.touchonly:
            return True
        else:
            for p in self.output_test:
                x = file.isnewer(self.input, p)
                if x:
                    logging.debug(f"Need to run because of file: {p}")
                    return x
            return x

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
            return True
        except subprocess.CalledProcessError as e:
            print("an error occured while executing {}".format(self.program))
            print(e.output)
            return False

    def checkOutput(self, file, comment="#"):
        """
        generic function to see if we have at least one none comment line
        in the output
        """
        with open(file) as f:
            for line in f:
                if not line.startswith(comment):
                    return True
        return False

    def cleanup(self, folders=None):
        """Remove temporary folders

        Given a list of folder names we will remove one after
        another, so that temporary directories are gone.
        """

        if folders is not None:
            if type(folders) == str:
                folders = [folders]
            for path in folders:
                if os.path.exists(path):
                    shutil.rmtree(path)
                    logging.debug("Removed temp folder: %s", path)
                else:
                    logging.debug("Cant remove folder, as it does not exists: %s", path)


# defining classes based on run for executing gmes and hmmscan
class gmes(run):
    def find_license(self):
        home = os.path.expanduser("~")
        pwd = os.getcwd()
        dirs = [home, pwd]
        key = ".gm_key"
        foundKey = False
        for directory in dirs:
            if os.path.exists(os.path.join(directory, key)):
                foundKey = True
                break
        return foundKey

    def run(self, cores=1):
        # before starting, check if the license key is anywhere we can find it
        if not self.find_license():
            logging.warning(
                "I could not find the license of GeneMark-ES. Make sure you have a valid license in your home directory with the name '.gm_key'. Will try to run GeneMark-ES regardless, but it will likely fail"
            )

        lst = [self.program, "-i", self.input, "-o", self.output, "-n", str(cores)]
        if self.touchonly:
            self.touch()
            return True
        try:
            subprocess.run(lst, check=True, shell=False)
            return True
        except subprocess.CalledProcessError:
            print("an error occured while executing {}".format(self.program))
            return False


class hmmpress(run):
    def run(self):
        # if sometimes we just touch, for debugging
        if self.touchonly:
            self.touch()
            return True
        lst = [self.program, "-f", self.input]
        try:
            subprocess.run(lst, check=True, shell=False, stdout=subprocess.DEVNULL)
            return True
        except subprocess.CalledProcessError:
            logging.error("an error occured while executing {}".format(self.program))
            exit(14)


class hmmer(run):
    def run(self, stdoutfile, hmmfiles, cores=1, modus="bitscore", evalue=False, training=False):
        # if sometimes we just touch, for debugging
        if self.touchonly:
            self.touch()
            return True

        # in training mode we dont cut of using ga but use an evalue
        # cutoff, so we can learn the bitscore tresholds
        if (training or modus == "evalue") and evalue:
            logging.debug("Using evlaue hmmsearch")
            lst = [
                self.program,
                "--cpu",
                str(cores),
                "-o",
                stdoutfile,
                "--tblout",
                self.output,
                "-E",
                str(evalue),
                hmmfiles,
                self.input,
            ]
        else:
            logging.debug("Using bitscore hmmsearch")
            lst = [
                self.program,
                "--cpu",
                str(cores),
                "-o",
                stdoutfile,
                "--tblout",
                self.output,
                "--cut_ga",
                hmmfiles,
                self.input,
            ]

        if self.debug:
            print(" ".join(lst))
        try:
            subprocess.run(lst, check=True, shell=False)
            return True
        except subprocess.CalledProcessError:
            logging.error("an error occured while executing {}".format(self.program))
            exit(16)

    def clean(self, hmmer, bedfile, resultfile, mindist=2000):
        bed = base.readbed(bedfile)
        fields = {"subject": 0, "evalue": 4, "score": 5, "profile": 2}
        table = []
        # intcoluns =
        ints = ["start", "stop"]
        # read in Hmmer results
        with open(hmmer) as h:
            for line in h:
                if line.startswith("#"):
                    continue
                l = line.split()
                n = {}
                for k, v in fields.items():
                    n[k] = l[v]
                for k in ["chrom", "start", "stop"]:
                    if n["subject"] in bed.keys():
                        n[k] = bed[n["subject"]][k]
                    else:
                        print("Could not find entry in bed file for {}".format(n["subject"]))
                # fix profile name
                pft = n["profile"].split(".")
                # if this profile is a subfamily, we want to name it FAMILY:Subfamily
                if pft[1].startswith("SF"):
                    n["profile"] = "{}:{}".format(pft[0], pft[1])
                else:
                    # else we just keep the family section
                    n["profile"] = pft[0]

                n["keep"] = True

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
            """
            returns true if row is within or close to the lastrow
            """
            # start within
            if row["start"] >= lastrow["start"] and row["start"] <= lastrow["stop"]:
                return True
            # end within
            if row["stop"] >= lastrow["start"] and row["stop"] <= lastrow["stop"]:
                return True

            # distance between any value with any other value:
            a = [row["start"], row["stop"]]
            b = [lastrow["start"], lastrow["stop"]]
            for v in a:
                for w in b:
                    if abs(v - w) <= distanceTresh:
                        return True

            return False

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
            if row["profile"] == lastrow["profile"] and row["chrom"] == lastrow["chrom"]:
                # check if both are close
                if closeOrWithin(row, lastrow, mindist):
                    if float(lastrow["evalue"]) < float(row["evalue"]):
                        keepI = False
                    else:
                        keepJ = False

            # save which row we keep
            if "keep" not in table[j].keys() or not keepJ:
                table[j]["keep"] = keepJ
            if "keep" not in table[i].keys() or not keepI:
                table[i]["keep"] = keepI

            # update j, so index to lastrow
            if not keepJ:
                j = i
            elif keepJ and keepI:
                j = i
            elif keepJ and keepI is False:
                # j stays the same
                j = j  # obv redundant, but nice for making sense of this code

        # write result to file
        cols = ["profile", "subject", "chrom", "start", "stop", "evalue", "score"]
        with open(resultfile, "w") as f:
            f.write("\t".join(cols + ["\n"]))
            for row in table:
                if not row["keep"]:
                    continue
                # print if we keep the row
                v = [str(row[k]) for k in cols] + ["\n"]
                f.write("\t".join(v))

        return resultfile


class hmmalign(run):
    def run(self, alignmentpath, hmmpath):
        # if sometimes we just touch, for debugging
        if self.touchonly:
            self.touch()
            return True

        lst = [
            self.program,
            "--outformat",
            "afa",
            "--mapali",
            alignmentpath,
            "--amino",
            hmmpath,
            self.input,
        ]
        # print(" ".join(lst))
        try:
            with open(self.output, "w") as f:
                subprocess.run(lst, check=True, shell=False, stdout=f)
            return True
        except subprocess.CalledProcessError:
            print("an error occured while executing {}".format(self.program))
            return False


class pplacer(run):
    def run(self, pkg, logfile, cores=1):
        # if sometimes we just touch, for debugging
        if self.touchonly:
            self.touch()
            return True

        lst = [
            self.program,
            "-o",
            self.output,
            "-p",
            "--keep-at-most",
            "5",
            "-m",
            "LG",
            "-j",
            str(cores),
            "-c",
            pkg,
            self.input,
        ]

        try:
            with open(logfile, "w") as lg:
                subprocess.run(lst, check=True, shell=False, stdout=lg)
            return True
        except subprocess.CalledProcessError:
            print("an error occured while executing {}".format(self.program))
            return False
        except Exception as e:
            print("Pplacer failed with unkown reason")
            print(e)
            return False

    def prepareAlignment(
        self, pplaceAlinment, hmmerOutput, proteinList, proteinFasta, config, cfg, tmpDir,
    ):
        # for each profile that we found a SCMG
        # make concatenated fasta (including new sequence)
        # align that to the hmm profile
        profiles = []
        with open(proteinList) as pl:
            for l in pl:
                profiles.append(l.strip())
        # read in hmmer output to get seqnames of target proteins
        # cols = []
        # hmmer = []
        # scmgsr = []

        # with open(hmmerOutput) as ho:
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
        scmgsr = [hit["profile"] for hit in hmmer]
        # scmgs = [p for p in set(scmgsr) if scmgsr.count(p) == 1]
        # db = [p for p in set(scmgsr) if scmgsr.count(p) == 2]
        scmgs = scmgsr
        # scmgs[-1] = db[0]
        scmgs.sort()
        # count the placements
        self.lenscmgs = len(scmgs)

        if self.lenscmgs > 0:
            # get the protein names for each scmgs
            proteinnames = []
            # for all profiles extract proteins
            hmmr = {}
            # load into dict
            for row in hmmer:
                hmmr[row["profile"]] = row["subject"]
            # make list of protein names
            for m in scmgs:
                proteinnames.append(hmmr[m])

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
            self.input = base.horizontalConcat(
                pplaceAlinment, alignments, scmgs, config.pkgfile("{}.refpkg".format("concat"), "aln_fasta"),
            )

    def reduceJplace(self, jplace, jplaceselection, placementCutoff=0.5):
        """
        reduces a jplace file to placements with at least placementCutoff
        post_prob
        """
        with open(jplace) as json_file:
            j = json.load(json_file)
            fields = j["fields"]
            newplacements = []
            for placement in j["placements"]:
                ps = placement["p"]
                kp = []
                for p in ps:
                    d = {k: v for k, v in zip(fields, p)}
                    if d["post_prob"] >= placementCutoff:
                        kp.append(p)

                if len(kp) > 0:
                    placementn = placement.copy()
                    placementn["p"] = kp
                    newplacements.append(placementn)
            jn = j.copy()
            jn["placements"] = newplacements

        with open(jplaceselection, "w") as f:
            json.dump(jn, f)
        return jplaceselection


class tog(run):
    def run(self):
        # if sometimes we just touch, for debugging
        logging.debug("Launching guppy tog")
        if self.touchonly:
            self.touch()
            return True

        lst = [self.program, "tog", self.input, "-o", self.output]
        try:
            subprocess.run(lst, check=True, shell=False)
            return True
        except subprocess.CalledProcessError:
            logging.error("an error occured while executing {}".format(self.program))
            return False


def checkDependencies(dep):
    """
    simple function to check all dependencies
    """
    allgood = True
    for name, program in dep.items():
        if run.which(program) is None:
            print("Dependency: Could not find {} ({})".format(program, name))
            allgood = False

    return allgood
