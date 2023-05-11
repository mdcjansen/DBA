#!/usr/bin/python3

import os, sys, logging, multiprocessing, re, glob, argparse, time, subprocess, threading, shutil
from datetime import datetime
#from Bio import SeqIO
#from Bio.SeqIO.QualityIO import FastqGeneralIterator

# Credentials
__author__ = "M.D.C. Jansen"
__version__ = "0.2"
__date__ = "28/04/2023"


# Check thread input
def cpu_threads(thread_input):
    if multiprocessing.cpu_count() > thread_input:
        return thread_input
    else:
        return multiprocessing.cpu_count()


# Validate existence of input
def valid_in(inputfolder):
    input = glob.glob(inputfolder)
    input = ''.join(inputfolder)
    if len(inputfolder) == 0:
        print("No input file found. Exiting programme.\n")
        logging.warning("No input file found. Exiting programme.")
        sys.exit(1)
    elif len(inputfolder) != 0:
        if os.path.exists(inputfolder) == True:
            print("\nFound folder: {fn}".format(fn= foldername))
            logging.info("Found folder {fn}".format(fn= foldername))
        elif os.path.exists(inputfolder) == False:
            print("\nNo input file found")
            logging.error("No input file found.")
            sys.exit(1)


# Checking output directory
def valid_out(outdir):
    if os.path.exists(outdir) == True:
        print("\nOutput directory already exists, do you want to clear the contents of the directory?")
        logging.warning("Output directory already exists, do you want to clear the contents of the directory?")
        answer_clearing = input("[y/n/exit]\n")
        if answer_clearing == "y" or answer_clearing == "Y":
            logging.info("Answer: {ac}".format(ac= answer_clearing))
            print("\nClearing directory...", end="\r")
            logging.info("Clearing directory...")
            os.chdir(outdir)
            os.system("find . ! -name 'Barcoding.log' ! -name '.' -type d -exec rm -rf {} +")
            os.system("find . ! -name 'Barcoding.log' ! -name '.' ! -type d -exec rm -rf {} +")
            print("Output directory has been cleared.\n")
            logging.info("Output directory has been cleared")
        elif answer_clearing == "n" or answer_clearing == "N":
            logging.info("Answer: {ac}".format(ac=answer_clearing))
            print("\nContinue analysis in selected folder or create default output directory for analysis?")
            logging.info("Continue analysis in selected folder or create default output directory for analysis?")
            answer_outdir = input("[con/def/exit]\n")
            if answer_outdir == "con" or answer_outdir == "CON":
                logging.info("Answer: {ao}".format(ao=answer_outdir))
                print("\nDirectory will not be cleared. Analyses resuming.\nNOTE: Existing data might be overwritten!\n")
                logging.warning("Directory will not be cleared. Analyses resuming. NOTE: Existing data might be overwritten!")
            elif answer_clearing == "def" or answer_clearing == "DEF":
                logging.info("Answer: {ao}".format(ao=answer_outdir))
                print("\nContinuing analysis with default output directory")
                logging.info("Continuing analysis with default output directory")
            elif answer_clearing == "exit" or answer_clearing == "EXIT":
                logging.info("Answer: {ac}".format(ac=answer_clearing))
                print("\nExiting analysis")
                logging.info("Analysis terminated by user, Analysis will be stopped")
                sys.exit(1)
            else:
                print("\nUnknown input. Please provide a valid input (con/CON - def/DEF - exit/EXIT). Terminating analysis.\n")
                logging.error("Unknown input. Please provide a valid input (con/CON - def/DEF - exit/EXIT). Terminating analysis.\n\n")
                sys.exit(1)
        elif answer_clearing == "exit" or answer_clearing == "EXIT":
            logging.info("Answer: {ac}".format(ac=answer_clearing))
            print("\nExiting analysis")
            logging.info("Analysis terminated by user, Analysis will be stopped")
            sys.exit(1)
        else:
            print("\nUnknown input. Please provide a valid input (y/Y - n/N - exit/EXIT). Terminating analysis.\n")
            logging.error("Unknown input. Please provide a valid input (y/Y - n/N - exit/EXIT). Terminating analysis.\n\n")
            sys.exit(1)
    elif os.path.exists(outdir) == False:
        os.makedirs(outdir)
        print("\nOutput directory: {od} has been created".format(od=outdir))
        logging.info("Output directory: {od} has been created".format(od=outdir))
        logging.info("Settings barcoding analysis:\n\nInput file:\t\t{fn}\nOutput directory:\t{od}\nThreads:\t\t{td}\nEstimated\nKeep all files:\t\t{kp}\n".format(fn=foldername, od=outdir, td=threads, kp=keep))

# Setup parser
def parser_config():
    global argument, parser, inputfolder, foldername, outdir, logfile, threads, keep, genbank, root, workdir
    parser = argparse.ArgumentParser(prog="Automated DNA barcoding analysis",
                                     description="ADba is designed to automate the process of DNA barcoding by utilising standard Sanger sequencing data and user provided reference gene numbers."
                                                 "This pipeline visualises sequences",
                                     usage="%(prog)s -i <inputfolder> [options]",
                                     epilog="")
    parser._optionals.title = "Arguments list"
    parser.add_argument("-v", "--version",
                        help="Prints program version",
                        action="version",
                        version="{v} {d} by {a}".format(v=__version__, d=__date__, a=__author__))
    parser.add_argument("-i",
                        metavar="[input]",
                        help="Input folder",
                        required=len(sys.argv) != 1)
    parser.add_argument("-g",
                        metavar="[genbank accession]",
                        help="genbank accession number for chromatogram annotation",
                        required= True)
    parser.add_argument("-o",
                        metavar="[output]",
                        help="Output directory",
                        required=False,
                        default="barcoding_output_{dt}/".format(dt=date.strftime("%d-%m-%Y_%H-%M-%S")))
    parser.add_argument("-t",
                        metavar="[threads]",
                        help="Amount of threads.",
                        required=False,
                        default=cpu_threads(8),
                        type=int)

    parser.add_argument("--keep",
                        metavar='',
                        help="Keep all files produced.",
                        required=False,
                        default=False,
                        type=bool,
                        nargs="?",
                        const=True)
    argument = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(0)
    root = os.path.dirname(os.path.realpath(__file__))
    inputfolder = os.path.abspath(argument.i)
    foldername = os.path.basename(argument.i)
    genbank = str(argument.g)
    outdir = os.path.abspath(argument.o)
    logfile = "{rt}/Barcoding.log".format(rt=root)
    threads = str(argument.t)
    keep = argument.keep
    workdir = os.path.join(root, "workdir")
    if not os.path.isdir(workdir):
        os.makedirs(workdir)
    if os.path.isdir(workdir):
        shutil.rmtree(workdir)
        os.makedirs(workdir)

# Fetching genbank database
def get_genbank():
    logging.info("Obtaining genbank entry")
    print("Obtaining genbank entry")
    genbankfunc = "bio fetch {gb} > {wd}/{gb}.gb".format(wd=workdir, gb=genbank)
    os.system(genbankfunc)
    logging.info("Acquired entry {gb} successfully".format(gb=genbank))
    print("Acquired entry {gb} successfully".format(gb=genbank))

# production of chromatogram and chromatogram score
def sangerseq():
    os.system("cp {inp}/*/*.ab1 {wd}/".format(inp=inputfolder, wd=workdir))
    sangerfunc = "sangerseq_viewer -s {wd}/{gb}.gb -q {wd}/ -o {o}/{cg}.png -l 200 -d 500"\
        .format(wd=workdir, gb=genbank, o=outdir, cg="chromatogram")
    os.system(sangerfunc)

def main():
    # Input validation and preparing log file
    parser_config()
    logging.basicConfig(#filename=logfile,
                        level=logging.DEBUG,
                        format="%(asctime)s - %(levelname)-8s - %(threadName)-10s - %(message)s",
                        handlers=[
                            logging.FileHandler(logfile),
                            logging.StreamHandler()
                            ]
                        )
    logging.info("Barcoding analysis initiated")
    valid_in(inputfolder)
    valid_out(outdir)

# RUNNING ANALYSIS ARGUMENTS

    get_genbank()
    sangerseq()

# Finishing analysis
    tend = int(time.time() - start_time)
    elapsed_time = "{:02d}:{:02d}:{:02d}".format(tend // 3600, (tend % 3600 // 60), tend % 60)
    bc_end = "\nAnalysis complete. The analysis took {time} to complete.".format(time=elapsed_time)
    print(bc_end)
    logging.info("Completed analysis\n{mg}\n\n\n".format(mg=bc_end))
    sys.exit(0)


if __name__ == '__main__':
    start_time = time.time()
    date = datetime.now()
    main()

