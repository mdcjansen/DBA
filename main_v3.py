#!/usr/bin/python3

import os, sys, logging, multiprocessing, re, glob, argparse, time, subprocess, threading, shutil
from datetime import datetime
#from Bio import SeqIO
#from Bio.SeqIO.QualityIO import FastqGeneralIterator

# Credentials
__author__ = "M.D.C. Jansen"
__version__ = "0.3"
__date__ = "08/05/2023"


# Check thread input
def cpu_threads(thread_input):
    if multiprocessing.cpu_count() > thread_input:
        return thread_input
    else:
        return multiprocessing.cpu_count()

# Validate input folder
def valid_in(inputfolder):
    if len(inputfolder) == 0:
        logging.warning("No input folder found. Exiting programme.")
        sys.exit(1)
    elif len(inputfolder) != 0:
        if os.path.exists(inputfolder) == True:
            logging.info("Found folder {fn}".format(fn=inputfolder))
        elif os.path.exists(inputfolder) == False:
            logging.error("No input folder found.")
            sys.exit(1)


# Validate output folder
def valid_out(outdir):
    if os.path.exists(outdir):
        logging.warning("Output directory already exists, do you want to clear the contents of the directory?")
        answer_clearing = input("[y/n/exit]: ")
        time.sleep(0.5)
        print("\033[A                             \033[A")
        if answer_clearing == "y" or answer_clearing == "Y":
            logging.info("Answer:\tYes")
            logging.info("Clearing directory...")
            os.chdir(outdir)
            logging.info("Output directory has been cleared")
        elif answer_clearing == "n" or answer_clearing == "N":
            logging.info("Answer:\tNo")
            logging.info("Continue analysis in selected folder or create default output directory for analysis?")
            answer_outdir = input("[continue/create/exit]: ")
            time.sleep(0.5)
            print("\033[A                             \033[A")
            if answer_outdir == "continue" or answer_outdir == "CONTINUE":
                logging.info("Answer:\tContinue")
                logging.warning("Directory will not be cleared. Analyses resuming. NOTE: Existing data might be overwritten!")
            elif answer_clearing == "create" or answer_clearing == "CREATE":
                logging.info("Answer:\tCreate")
                logging.info("Continuing analysis with default output directory")
            elif answer_clearing == "exit" or answer_clearing == "EXIT":
                logging.info("Answer:\tExit")
                tend = int(time.time() - start_time)
                elapsed_time = "{:02d}:{:02d}:{:02d}".format(tend // 3600, (tend % 3600 // 60), tend % 60)
                logging.info("Analysis terminated by user. Analysis stopped after: {et}\n\n\n".format(et=elapsed_time))
                sys.exit(1)
            else:
                logging.error("Unknown input. Please provide a valid input (con/CON - def/DEF - exit/EXIT). Terminating analysis.\n\n")
                sys.exit(1)
        elif answer_clearing == "exit" or answer_clearing == "EXIT":
            logging.info("Answer: {ac}".format(ac=answer_clearing))
            tend = int(time.time() - start_time)
            elapsed_time = "{:02d}:{:02d}:{:02d}".format(tend // 3600, (tend % 3600 // 60), tend % 60)
            logging.info("Analysis terminated by user. Analysis stopped after: {et}\n\n\n".format(et=elapsed_time))
            sys.exit(1)
        else:
            logging.error("Unknown input. Please provide a valid input (y/Y - n/N - exit/EXIT). Terminating analysis.\n\n")
            tend = int(time.time() - start_time)
            elapsed_time = "{:02d}:{:02d}:{:02d}".format(tend // 3600, (tend % 3600 // 60), tend % 60)
            logging.error("Analysis terminated due to bad input. Analysis stopped after: {et}\n\n\n".format(et=elapsed_time))
            sys.exit(1)
    elif not os.path.exists(outdir):
        os.makedirs(outdir)
        logging.info("Output directory: {od} has been created".format(od=outdir))


# Setup parser
def parser_config():
    global argument, parser, inputfolder, outdir, logfile, threads, keep, genbank, root, workdir
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
                        default=cpu_threads(20),
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
    genbank = str(argument.g)
    outdir = os.path.abspath(argument.o)
    logfile = "{rt}/barcoding.log".format(rt=root)
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
    cmd = "bio fetch {gb} > {wd}/{gb}.gb".format(wd=workdir, gb=genbank)
    proc = subprocess.run([cmd],
                          shell=True,
                          text=True,
                          capture_output=True
                          )
    if len(proc.stderr) != 0:
        logging.error(proc.stderr)
        tend = int(time.time() - start_time)
        elapsed_time = "{:02d}:{:02d}:{:02d}".format(tend // 3600, (tend % 3600 // 60), tend % 60)
        logging.error("Analysis terminated after: {et}\n\n\n".format(et=elapsed_time))
        sys.exit(1)
    else:
        logging.info("Genbank entry: {gb} obtained successfully".format(gb=genbank))


# production of chromatogram
def sangerseq():
    logging.info("Creating chromatogram")
    os.system("cp {inp}/*/*.ab1 {wd}/".format(inp=inputfolder, wd=workdir))
    cmd= "sangerseq_viewer -s {wd}/{gb}.gb -q {wd}/ -o {o}/{cg}.png -l 200 -d 100"\
        .format(wd=workdir, gb=genbank, o=outdir, cg="chromatogram")
    proc = subprocess.run([cmd],
                          shell=True,
                          text=True,
                          capture_output=True
                          )
    if len(proc.stderr) != 0:
        logging.error(proc.stderr)
        tend = int(time.time() - start_time)
        elapsed_time = "{:02d}:{:02d}:{:02d}".format(tend // 3600, (tend % 3600 // 60), tend % 60)
        logging.error("Analysis terminated after: {et}\n\n\n".format(et=elapsed_time))
        sys.exit(1)
    else:
        logging.info("Successfully produced chromatogram")
    # -d to 500

def main():
    # Input validation and preparing log file
    parser_config()
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(levelname)-8s - %(threadName)-10s - %(message)s",
                        handlers=[
                            logging.FileHandler(logfile),
                            logging.StreamHandler()
                            ]
                        )
    logging.info("Barcoding analysis initiated")
    valid_in(inputfolder)
    valid_out(outdir)
    logging.info("Settings barcoding analysis:\n\n"
                 "Input folder:\t\t{inp}\n"
                 "Output directory:\t{od}\n"
                 "Genbank accession no.:\t{gb}\n"
                 "Threads:\t\t{td}\n"
                 "Keep all files:\t\t{kp}\n"
                 .format(inp=inputfolder, od=outdir, gb=genbank, td=threads, kp=keep))

# RUNNING ANALYSIS ARGUMENTS

    get_genbank()
    sangerseq()

# Finishing analysis
    tend = int(time.time() - start_time)
    elapsed_time = "{:02d}:{:02d}:{:02d}".format(tend // 3600, (tend % 3600 // 60), tend % 60)
    logging.info("Analysis successful. Time to completion: {et}".format(et=elapsed_time))
    sys.exit(0)


if __name__ == '__main__':
    start_time = time.time()
    date = datetime.now()
    main()

