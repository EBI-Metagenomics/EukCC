#
# file to keep track of simple function

import os
import datetime

def exists(f):
    return(os.path.exists(f))


def log(m, verbose = True):
    ts = datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S")
    if verbose:
        print("{}: {}".format(ts, m))
        
  