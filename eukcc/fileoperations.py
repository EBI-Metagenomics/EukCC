#!/usr/bin/env python
#
# provides all file operation functions
# used inthis package
#
import os
import logging


class file():

    @staticmethod
    def isfile(f):
        return(os.path.exists(f))

    @staticmethod
    def exists(f):
        return(os.path.exists(f))

    @staticmethod
    def isdir(d, create=True):
        """
        check if dir exists and create if not
        """
        if not os.path.isdir(d):
            if create:
                try:
                    os.makedirs(d)
                    return(True)
                except OSError as e:
                    logging.warning(f"Could not create dir: {d}\n{e}")
                    return(False)
            else:
                return(False)
        else:
            return(True)

    @staticmethod
    def isnewer(fileA, fileB):
        '''
        Check if file A is newer than file B
        '''
        if not file.isfile(fileA):
            logging.warning("{} is not a file".format(fileA))
            return(False)

        if not file.isfile(fileB):
            # return true bc fileB does not exists!
            return(True)
        return(os.stat(fileA).st_mtime > os.stat(fileB).st_mtime)

    @staticmethod
    def touch(path):
        with open(path, 'a'):
            os.utime(path, None)




