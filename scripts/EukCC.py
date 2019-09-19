#!/usr/bin/env python
import argparse
import os
import glob
from eukcc.eukcc import eukcc as runeukcc
from eukcc.base import log
from eukcc.base import log
from eukcc.fileoperations import file


# set arguments
# arguments are passed to classes
parser = argparse.ArgumentParser(description='Evaluate completeness and contamination of a MAG')
parser.add_argument('outdir', type=str, default="./",
                    help="Location for the output. Names will be prefixed using \
                          the bin filenames")
parser.add_argument('fasta',  type=str, nargs='+',
                    help='run script on this bin (fasta file)')
parser.add_argument('--configdir','-c', type=str, metavar="PATH",
                    default="./",
                    help='directory containing set sets and hmms')
parser.add_argument('--ncores','-n', metavar="int", type=int,
                    default=1,
                    help='set number of cores for GeneMark-ES, pplacer and Hmmer')
parser.add_argument('--hmm', dest='hmm',  type=str, 
                    default=None, help='run hmmer on all these HMMs instead')
parser.add_argument('--training', dest='training', action='store_true', 
                    default=False, help='run EukCC in training mode (needed to create a new release of the DB)')
parser.add_argument('--evalue', dest='evalue',  type=float, 
                    default=1e-5, help='use this evalue cutoff for hmmer')
parser.add_argument('--bed','-b', metavar="file.bed", type=str,
                    default=None,
                    help='pass bedfile if you called genes manually. \
                    Assumes only a single fasta (protein) is passed and implies --noglob')
parser.add_argument('--force', '-f', dest='force', action='store_true',
                    default=False, help='force rerun of computation even if \
                                          output is newer than input. Don\'t resume previous run.')
parser.add_argument('--fplace', '-p', dest='fplace', action='store_true',
                    default=False, help='force rerun of placement and subsequent steps')
parser.add_argument('--noplace', dest='noplace', action='store_true',
                    default=False, help='Do not place MAG, stop after running GeneMark-ES')
parser.add_argument('--noglob', '-g', dest='noglob', action='store_true',
                    default=False, help='Do not expand paths using glob')
parser.add_argument('--quiet', '-q', dest='quiet', action='store_true',
                    default=False, help='silcence most output')
parser.add_argument('--debug', '-d',  action='store_true',
                    default=False, help='debug and thus ignore safety')
args = parser.parse_args()


###############################################
# starting the analysis
log("Running eukcc for {} bin{}".format(len(args.fasta), "s" if len(args.fasta) > 1 else ""))

# create output if not exists
if not file.isdir(args.outdir):
    exit()

# check if a protein fasta was passed (implied )
if args.bed is not None:
    # set no glob
    args.noglob = True
    args.isprotein = True
else:
    args.isprotein = False
    
    
# check if we can expand glob:
if len(args.fasta) == 1 and not args.noglob:
    log("Expanding paths using glob", not args.quiet)
    args.fasta = glob.glob(args.fasta[0])

for fa in args.fasta:
    # check if input file exists if not break
    if not file.isfile(fa):
        log("Error: Could not find fasta:\n{}".format(fa))
        continue
    
    # get fasta name
    name = (os.path.splitext(os.path.basename(fa))[0])
    
    if len(args.fasta) > 1:
        print("")
        log("Running EukCC for {}".format(name))
    
    
    try:
        runeukcc(fa, args.configdir, 
                 outdir = os.path.join(args.outdir, name),
                 threads = args.ncores,
                 evalue = args.evalue,
                 force = args.force,
                 fplace = args.fplace,
                 isprotein = args.isprotein,
                 bedfile = args.bed,
                 hmm = args.hmm,
                 debug = args.debug,
                 noplace = args.noplace,
                 training = args.training)
    except Exception as e:
        log("Could not run EukCC for {}\n check logs for details".format(name))
        print(e)
        
    if len(args.fasta) > 1:
        print("")
        


###############################################
