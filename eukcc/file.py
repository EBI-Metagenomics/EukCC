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
#
import os
import logging
import shutil


class file:
    @staticmethod
    def isfile(f):
        return os.path.exists(f)

    @staticmethod
    def exists(f):
        return os.path.exists(f)

    @staticmethod
    def isdir(d):
        """
        check if dir exists and create if not
        """
        if not os.path.isdir(d):
            os.makedirs(d)
            return os.path.isdir(d)
        else:
            return os.path.isdir(d)

    @staticmethod
    def isnewer(fileA, fileB):
        """
        Check if file A is newer than file B
        """
        if not file.isfile(fileA):
            raise OSError("Could not find {}".format(fileA))

        if not file.isfile(fileB):
            # return true bc fileB does not exists!
            return 2
        return os.stat(fileA).st_mtime > os.stat(fileB).st_mtime

    @staticmethod
    def touch(path):
        with open(path, "a"):
            os.utime(path, None)

    @staticmethod
    def remove(path):
        if os.path.exists(path):
            if os.path.isdir(path):
                shutil.rmtree(path)
            else:
                os.remove(path)

    @staticmethod
    def delete_but(folder, keep=None):
        if type(keep) is str:
            keep = [keep]
        if keep is None:
            keep = []
        # list files
        fls = os.listdir(folder)
        for f in fls:
            if f not in keep:
                file.remove(os.path.join(folder, f))
