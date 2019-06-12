#!/usr/bin/env python
import datetime

def log(m, quiet = False):
    ts = datetime.datetime.now().strftime("%d/%m/%Y, %H:%M:%S")
    if not quiet:
        print("{}: {}".format(ts, m))
