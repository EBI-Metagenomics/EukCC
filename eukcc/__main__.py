import logging
import configargparse
from eukcc import workflow
import os


    
def main():
    # set arguments
    # arguments are passed to classes
    parser = configargparse.ArgumentParser(description='Evaluate completeness \
                        and contamination of a MAG.')
    parser.add_argument('fasta',  type=str,
                        help='run script on this bin (fasta file)')
    parser.add_argument('--db', type=str, required=True,
                        help='Path to EukCC DB')
    parser.add_argument('--outdir', '-o', type=str, default="./",
                        help="Location for the output. Names will be prefixed using \
                              the bin filenames")
    parser.add_argument('--config', '-c', type=str,
                        required=False, is_config_file=True,
                        help='Config file to define parameters, YAML')
    parser.add_argument('--ncores', '-n', metavar="int", type=int,
                        default=1,
                        help='set number of cores for GeneMark-ES, pplacer and Hmmer')
    parser.add_argument('--hmm', dest='hmm',  type=str, 
                        default=None, help='run hmmer on all these HMMs instead')
    parser.add_argument('--training', dest='training', action='store_true', 
                        default=False, help='run EukCC in training mode (needed to create a new release of the DB)')
    parser.add_argument('--bed', '-b', metavar="file.bed", type=str,
                        default=None,
                        help='pass bedfile if you called genes manually. \
                        Assumes only a single fasta (protein) is passed and implies --noglob')
    parser.add_argument('--force', '-f', dest='force', action='store_true',
                        default=False, help='force rerun of computation even if \
                                              output is newer than input. Don\'t resume previous run.')
    parser.add_argument('--fplace', '-p', dest='fplace', action='store_true',
                        default=False, help='force rerun of placement and subsequent steps')
    parser.add_argument('--noglob', '-g', dest='noglob', action='store_true',
                        default=False, help='Do not expand paths using glob')
    parser.add_argument('--quiet', '-q', dest='quiet', action='store_true',
                        default=False, help='silcence most output')
    parser.add_argument('--debug', '-d',  action='store_true',
                        default=False, help='debug and thus ignore safety')
    parser.add_argument('--HPA',  default=False, action='store_true',
                        help = "Set placement method to HPA")
    parser.add_argument('--nPlacements', type=int, default=2, metavar = "n",
                        help = "Set number of proteins to support location \
                                in tree (default: 2)")
    parser.add_argument('--fulllineage', default = False, action='store_true',
                        help = "Output full lineage for MAGs")
    parser.add_argument('--minPlacementLikelyhood', default = 0.4, type = float,
                        metavar="float",
                        help = "minimal pplacer likelyhood (default: 0.4)")
    parser.add_argument('--mindist', type=int, default=2000, metavar="n",
                        help = "Distance to collapse hits (default: 2000)")
    parser.add_argument('--touch', default=False, action='store_true',
                        help="Do not run, but touch all output files")
    options = parser.parse_args()

    # define logging
    logLevel = logging.INFO
    if options.quiet:
        logLevel = logging.WARNING
    elif options.debug:
        logLevel = logging.DEBUG
    logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %H:%M:%S: ',
                        level=logLevel)
    logging.debug("Launching EukCC in debug mode")
    logging.info("Starting EukCC")

    # Now we start the run with EukCC
    # All magic numbers should be defined in info.py if they are not
    # part of the configuration options
    m = workflow.eukcc(options)
    
    # skip gene predition if this is already protein sequences
    if options.bed is None and not m.stopnow():
        # run gmes
        proteinfaa, bedfile = m.gmes(options.fasta)
    else:
        proteinfaa = options.fasta
        bedfile = options.bed
    
    # run hmm file if we are asked to
    # this is needed during for training 
    if m.cfg['training'] or m.cfg['hmm']:
        logging.info("Running on custom hmm for training mode")
        m.runPlacedHMM(m.cfg['hmm'], proteinfaa, bedfile)
        logging.info("Stopping now as we are only doing training")
        exit()
    
    # place using pplacer and hmmer
    m.placed = m.place(proteinfaa, bedfile)
    
    # concat hmms for hmmer
    hmmfile = m.concatHMM(m.placed)
    # run Hmmer for sets of placement
    hits = m.runPlacedHMM(hmmfile, proteinfaa, bedfile)
    # infer lineage
    _ = m.inferLineage(m.placed[m.cfg['placementMethod']])
    _ = m.plot()
        
    # estimate completeness and contamiantion
    outputfile = os.path.join(m.cfg['outdir'], m.cfg['outfile'])
    m.estimate(hits, outputfile, m.placed[m.cfg['placementMethod']])
