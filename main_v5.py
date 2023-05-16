#!/usr/bin/python3

import argparse
import logging
import multiprocessing
import os
import shutil
import subprocess
import sys
import time
from datetime import datetime

# Credentials
__author__ = "M.D.C. Jansen"
__version__ = "0.5.5"
__date__ = "16/05/2023"


# Configure inputs
def config():
    global argument, parser, inputfolder, outdir, logfile, threads, keep, genbank, root, workdir

    # Setup parser
    parser = argparse.ArgumentParser(prog="Automated DNA barcoding analysis",
                                     description="ADba is designed to automate the process of DNA barcoding "
                                                 "by utilising standard Sanger sequencing data "
                                                 "and user provided reference gene numbers."
                                                 "This pipeline visualises sequences OTHER THINGS IT DOES",
                                     usage="%(prog)s -i <inputfolder> [options]",
                                     epilog="")
    parser._optionals.title = "Arguments list"
    parser.add_argument("-v", "--version",
                        help="Prints program version",
                        action="version",
                        version="Version: {v} Date: {d} By {a}".format(v=__version__, d=__date__, a=__author__))
    parser.add_argument("-i",
                        metavar="[input]",
                        help="Input folder",
                        required=len(sys.argv) != 1)
    parser.add_argument("-g",
                        metavar="[genbank accession]",
                        help="genbank accession number for chromatogram annotation",
                        required=len(sys.argv) != 1)
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
    genbank = str(argument.g)
    threads = str(argument.t)
    keep = argument.keep

    # Configuring logfile
    logfile = "{rt}/barcoding.log".format(rt=root)
    if os.path.isfile(logfile):
        os.remove(logfile)
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(levelname)-8s - %(message)s",
                        handlers=[
                            logging.FileHandler(logfile),
                            logging.StreamHandler()
                            ]
                        )
    logging.info("Barcoding analysis initiated")

    # Validating input
    inputfolder = os.path.abspath(argument.i)
    if len(inputfolder) == 0:
        logging.warning("No input folder found.")
        tend = int(time.time() - start_time)
        elapsed_time = "{:02d}:{:02d}:{:02d}".format(tend // 3600, (tend % 3600 // 60), tend % 60)
        logging.error("Analysis terminated after: {et}\n\n\n".format(et=elapsed_time))
        sys.exit(1)
    elif len(inputfolder) != 0:
        if os.path.exists(inputfolder):
            logging.info("Found folder {fn}".format(fn=inputfolder))
            global folder_count, folder_names
            folder_names = os.listdir(inputfolder)
            folder_count = len([dir for dir in os.listdir(inputfolder) if os.path.isdir(os.path.join(inputfolder, dir))])
        elif not os.path.exists(inputfolder):
            logging.error("No input folder found.")
            tend = int(time.time() - start_time)
            elapsed_time = "{:02d}:{:02d}:{:02d}".format(tend // 3600, (tend % 3600 // 60), tend % 60)
            logging.error("Analysis terminated after: {et}\n\n\n".format(et=elapsed_time))
            sys.exit(1)

    # Validating output
    outdir = os.path.abspath(argument.o)
    if os.path.exists(outdir):
        logging.warning("Output directory already exists, Continuing analysis in default output folder")
        outdir = "barcoding_output_{dt}/".format(dt=date.strftime("%d-%m-%Y_%H-%M-%S"))
        os.makedirs(outdir)
        logging.info("Output directory: {od} has been created".format(od=outdir))
    elif not os.path.exists(outdir):
        os.makedirs(outdir)
        logging.info("Output directory: {od} has been created".format(od=outdir))

    # Creating working directory
    workdir = os.path.join(root, "workdir")
    if not os.path.isdir(workdir):
        os.makedirs(workdir)
    if os.path.isdir(workdir):
        shutil.rmtree(workdir)
        os.makedirs(workdir)


# Check thread input
def cpu_threads(thread_input):
    if multiprocessing.cpu_count() > thread_input:
        return thread_input
    else:
        return multiprocessing.cpu_count()


# Running analysis
def process_run(cmd_in, process_start, completion):
    logging.info(process_start)
    process = subprocess.run([cmd_in],
                             shell=True,
                             text=True,
                             capture_output=True
                             )
    if len(process.stderr) != 0:
        if cmd_in == cmd_muscle or cmd_in == cmd_muscle_data:
            pass
        else:
            process_err(process)
    else:
        logging.info(completion)


# Error posting during analysis
def process_err(process):
    logging.error(process.stderr)
    tend = int(time.time() - start_time)
    elapsed_time = "{:02d}:{:02d}:{:02d}".format(tend // 3600, (tend % 3600 // 60), tend % 60)
    logging.error("Analysis terminated after: {et}\n\n\n".format(et=elapsed_time))
    sys.exit(1)


# Preparing working folder for processing of  datasets
def work_prep():
    logging.info("Preparing working folder for {fc} inputs"
                 .format(fc=str(folder_count)))
    num_folder = 0
    while num_folder != folder_count:
        os.makedirs("{wd}/{fn}"
                    .format(wd=workdir, fn=folder_names[num_folder]))
        os.system("cp {inp}/{fn}/*/* {wd}/{fn}/"
                  .format(inp=inputfolder, fn=folder_names[num_folder], wd=workdir))
        num_folder += 1


# Multi-threading analysis of single-threaded chromatogram analysis
def multi_sanger():
    processes = []
    for i in range(0, len(folder_names)):
        cmd = "sangerseq_viewer -s {wd}/{gb}.gb -q {wd}/{fn} -o {o}/{cg}_{fn}.png -l 200 -d 100" \
            .format(wd=workdir, gb=genbank, fn=folder_names[int(i)], o=outdir, cg="chromatogram")
        st_log = "Creating chromatogram for {fn}".format(fn=folder_names[int(i)])
        ed_log = "Successfully produced chromatogram for {fn}".format(fn=folder_names[int(i)])
        p = multiprocessing.Process(target=process_run, args=(cmd, st_log, ed_log))
        processes.append(p)
        p.start()
    for j in processes:
        j.join()


def main():
    # Parser setup, input/output validation, and preparing log file
    global cmd_muscle, cmd_muscle_data
    config()
    logging.info("Settings barcoding analysis:\n\n"
                 "Input folder:\t\t{inp}\n"
                 "Output directory:\t{od}\n"
                 "Genbank accession no.:\t{gb}\n"
                 "Threads:\t\t{td}\n"
                 "Keep all files:\t\t{kp}\n"
                 .format(inp=inputfolder, od=outdir, gb=genbank, td=threads, kp=keep))

    # Analysis inputs
    cmd_genbank = "bio fetch {gb} > {wd}/{gb}.gb ; " \
                  "bio fetch {gb} --format fasta > {wd}/{gb}.fasta"\
        .format(wd=workdir, gb=genbank)
    st_genbank = "Obtaining genbank entry: {gb}"\
        .format(gb=genbank)
    ed_genbank = "Successfully obtained Genbank entry: {gb}"\
        .format(gb=genbank)
    cmd_muscle = "cat {wd}/*/*.fa > {wd}/{sq} && " \
                 "muscle -align {wd}/{sq} -output {od}/{ou} -threads {ts} && " \
                 "muscle -align {wd}/{sq} -output {od}/{ot} -threads {ts} -diversified"\
        .format(wd=workdir, sq="combined.fa", od=outdir, ou="aln.fa", ot="diversified_aln_confseq.efa", ts=threads)
    st_muscle = "Starting alignment with muscle"
    ed_muscle = "Successfully performed alignment"
    cmd_muscle_data = "{me} -addconfseq {od}/{it} -output {od}/{ou} ; " \
                      "{me} -letterconf {od}/{it} -ref {od}/{al} -output {od}/{ot}" \
                      " -html {od}/{hl} -jalview {od}/{jv} ; " \
                      "{me} -efatats {od}/{it} -log {od}/efastats.log ; " \
                      "{me} -disperse {od}/{it} -log {od}/disperse.log "\
        .format(me="muscle", od=outdir, it="/diversified_aln.efa", ou="diversified_aln_confseq.efa",
                ot="letterconf.afa", al="aln.fa", hl="letterconf.html", jv="letterconf_jalview.features")
    st_muscle_data = "Starting collection of additional muscle alignment data"
    ed_muscle_data = "Successfully obtained additional data"

    # Running analysis
    process_run(cmd_genbank, st_genbank, ed_genbank)
    work_prep()
    multi_sanger()
    process_run(cmd_muscle, st_muscle, ed_muscle)
    process_run(cmd_muscle_data, st_muscle_data, ed_muscle_data)

    # Finishing analysis
    tend = int(time.time() - start_time)
    elapsed_time = "{:02d}:{:02d}:{:02d}".format(tend // 3600, (tend % 3600 // 60), tend % 60)
    logging.info("Analysis successful. Time to completion: {et}".format(et=elapsed_time))
    sys.exit(0)


if __name__ == '__main__':
    start_time = time.time()
    date = datetime.now()
    main()
