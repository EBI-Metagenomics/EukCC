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
import json
import logging
import yaml
from eukcc import base
from eukcc.base import log


# default magic numbers are saved here
# please keep it sorted alphabetically and coument
defaults = {
    "evalue": 1e-5,  # not used evalue
    "nEvals": 3,  # per placement show top n inferrals
    "minProfiles": 20,  # minimal profiles in set
    "minSupport": 2,  # minimal profiles to support set
    "trainingEvalue": 10,  # evalue used in training mode
}


class eukinfo:
    def __init__(self, options):
        self.setConfig(options)

    def setConfig(self, options):
        self.cfg = defaults
        # loop over options and import into our configs
        for k, v in vars(options).items():
            self.cfg[k] = v

        self.loadDBinfo()

        # placements Methis
        if options.HPA:
            self.cfg["placementMethod"] = "HPA"
        else:
            self.cfg["placementMethod"] = "LCA"
        # figure if pplacer cores need to be adjusted
        if self.cfg["ncorespplacer"] < 1:
            logging.debug("Set pplacer cores to the same as all others (%d)", self.cfg["ncores"])
            self.cfg["ncorespplacer"] = self.cfg["ncores"]

        # define location of placement HMMs
        self.placementHMMs = os.path.join(options.db, "hmms/concat.hmm")
        self.tree = self.pkgfile("concat.refpkg", "tree")

    def loadDBinfo(self):
        dbinfopath = os.path.join(self.cfg["db"], "dbconfig.yml")
        if not os.path.exists(dbinfopath):
            logging.info("This is an old database, it does not profide machine readable information about the database")
            # for now we maintain this section, setting defaults
            # to make the software compatible with testing old DB versions
            # but in the future this should be changed
            self.cfg["dbinfo"] = {}
            self.cfg["dbinfo"]["modus"] = "bitscore"
            logging.debug(
                "Remove this manual setting of options to break backwards compatibility,\
                        to make use of old DB versions impossible"
            )
        else:
            with open(dbinfopath, "r") as stream:
                try:
                    dbinfo = yaml.safe_load(stream)
                    # copy the dbinfo into the config
                    self.cfg["dbinfo"] = {}
                    for key, value in dbinfo.items():
                        self.cfg["dbinfo"][key] = value
                except yaml.YAMLError as exc:
                    logging.error("We could not open the DBinfo file:")
                    print(exc)

    def checkForFiles(self, dirname):
        required = ["profile.list", "refpkg", "hmms/concat.hmm", "sets/setinfo.csv"]
        for f in required:
            p = os.path.join(dirname, f)
            if not base.exists(p):
                print("Configuartion folder does not contain: {}".format(f))
                return False
        return True

    def pkgfile(self, name, t):
        """
        get a file path for a refpkg package
        """
        info = self.readInfo(name)
        p = os.path.join(self.cfg["db"], "refpkg", name, info["files"][t])
        if base.exists(p):
            return p
        else:
            log("Could not find: {}".format(p))
            exit(12)

    def readInfo(self, name):
        p = os.path.join(self.cfg["db"], "refpkg", name, "CONTENTS.json")
        # raise error if we cant find the file
        if not base.exists(p):
            log("Could not find {}".format(p))
            exit(13)
        # read and return json
        with open(p) as json_file:
            j = json.load(json_file)
            return j
