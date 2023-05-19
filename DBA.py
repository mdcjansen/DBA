#!/usr/bin/python3

import argparse
import glob
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
__version__ = "0.6.1"
__date__ = "19/05/2023"


# Configure inputs
def config():
    global argument, parser, inputfolder, outdir, logfile, threads, keep, genbank_nucl, genbank_prot, root, workdir

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
                        metavar="[genbank nucleotide accession]",
                        help="genbank nucleotide accession number for chromatogram annotation",
                        required=len(sys.argv) != 1)
    parser.add_argument("-p",
                        metavar="[genbank protein accession]",
                        help="genbank protein accession number for chromatogram annotation",
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
    genbank_nucl = str(argument.g)
    genbank_prot = str(argument.p)
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
            logging.info("Found folder: {fn}".format(fn=inputfolder))
            global folder_count, folder_names, inputfasta
            folder_names = [dir for dir in os.listdir(inputfolder) if os.path.isdir(os.path.join(inputfolder, dir))]
            folder_count = len(folder_names)
            inputfiles = [file for file in os.listdir(inputfolder) if file.endswith(".fasta")
                          if os.path.isfile(os.path.join(inputfolder, file))]
            if len(inputfiles) == 1:
                inputfasta = inputfiles[0]
                logging.info("Found file: {inf}"
                             .format(inf=inputfasta))
            else:
                logging.error("Incorrect number of input FASTA files found."
                              "Expected one, found {nf}".format(nf=len(inputfiles)))
                logging.error("No input folder found.")
                tend = int(time.time() - start_time)
                elapsed_time = "{:02d}:{:02d}:{:02d}".format(tend // 3600, (tend % 3600 // 60), tend % 60)
                logging.error("Analysis terminated after: {et}\n\n\n".format(et=elapsed_time))
                sys.exit(1)
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

    # Creating and preparing working directory
    logging.info("Preparing working folder for {fc} inputs"
                 .format(fc=str(folder_count)))
    workdir = os.path.join(root, "workdir")
    if not os.path.isdir(workdir):
        os.makedirs(workdir)
    if os.path.isdir(workdir):
        shutil.rmtree(workdir)
        os.makedirs(workdir)
    num_folder = 0
    while num_folder != folder_count:
        os.makedirs("{wd}/{fn}"
                    .format(wd=workdir, fn=folder_names[num_folder]))
        os.system("cp {inp}/{fn}/*/* {wd}/{fn}/"
                  .format(inp=inputfolder, fn=folder_names[num_folder], wd=workdir))
        num_folder += 1
    shutil.copy(inputfasta, workdir)


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
    if process.stderr is not None and len(process.stderr) != 0:
        logging.error(process.stderr)
        tend = int(time.time() - start_time)
        elapsed_time = "{:02d}:{:02d}:{:02d}".format(tend // 3600, (tend % 3600 // 60), tend % 60)
        logging.error("Analysis terminated after: {et}\n\n\n".format(et=elapsed_time))
        sys.exit(1)
    else:
        logging.info(completion)


# Multi-threading analysis of single-threaded chromatogram analysis
def multi_sanger():
    processes = []
    for i in range(0, len(folder_names)):
        cmd = "sangerseq_viewer -s {wd}/{gbn}.gb -q {wd}/{fn} -o {o}{cg}_{fn}.pdf" \
            .format(wd=workdir, gbn=genbank_nucl, fn=folder_names[int(i)], o=outdir, cg="chromatogram")
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
                 "Input folder:\t\t\t\t{inp}\n"
                 "Output directory:\t\t\t{od}\n"
                 "Genbank nucleotide accession no.:\t{gbn}\n"
                 "Genbank protein accession no.:\t\t{gbp}\n"
                 "Threads:\t\t\t\t{td}\n"
                 "Keep all files:\t\t\t\t{kp}\n"
                 .format(inp=inputfolder, od=outdir, gbn=genbank_nucl, gbp= genbank_prot, td=threads, kp=keep))

    # Analysis inputs

    # Obtaining data from genbank
    cmd_genbank = "bio fetch {gbn} > {wd}/{gbn}.gb ; " \
                  "bio fetch {gbn} --format fasta > {wd}/{gbn}.fasta ; " \
                  "bio fetch {gbp} --format fasta > {wd}/{gbp}.prot"\
        .format(wd=workdir, gbn=genbank_nucl, gbp=genbank_prot)
    st_genbank = "Obtaining genbank entries:\t\t\t{gbn},\t{gbp}"\
        .format(gbn=genbank_nucl, gbp=genbank_prot)
    ed_genbank = "Successfully obtained Genbank entries:\t{gbn},\t{gbp}"\
        .format(gbn=genbank_nucl, gbp=genbank_prot)

    # Running BLASTn and BLASTx
    cmd_blast_db = "makeblastdb -in {wd}/{gbn}.fasta -dbtype nucl -out {wd}/{gbn}_nucldb ; " \
                   "makeblastdb -in {wd}/{gbp}.prot -dbtype prot -out {wd}/{gbp}_protdb"\
        .format(wd=workdir, gbn=genbank_nucl, gbp=genbank_prot)
    st_blast_db = "Creating blast databases"
    ed_blast_db = "Successfully created databases"
    cmd_blast_run = "blastn -query {wd}/{inf} -db {wd}/{gbn}_nucldb -out {ou}blastn.tsv -outfmt 6 ; " \
                    "blastx -query {wd}/{inf} -db {wd}/{gbp}_protdb -out {ou}blastx.tsv -outfmt 6"\
        .format(wd=workdir, inf=inputfasta, ou=outdir, gbn=genbank_nucl, gbp=genbank_prot)
    st_blast_run = "Starting BLASTn and BLASTx analysis"
    ed_blast_run = "Successfully completed BLASTn and BLASTx analysis"

    # Aligning sequences with MUSCLE
    cmd_muscle = "muscle -align {wd}/{sq} -output {od}{ou} -threads {ts} && " \
                 "muscle -align {wd}/{sq} -output {od}{ot} -threads {ts} -diversified"\
        .format(wd=workdir, sq="combined.fa", od=outdir, ou="aln.fa", ot="diversified_aln_confseq.efa", ts=threads)
    st_muscle = "Starting alignment with muscle"
    ed_muscle = "Successfully performed alignment"
    cmd_muscle_data = "{me} -addconfseq {od}/{it} -output {od}/{ou} ; " \
                      "{me} -letterconf {od}/{it} -ref {od}/{al} -output {od}/{ot}" \
                      " -html {od}{hl} -jalview {od}/{jv} ; " \
                      "{me} -efatats {od}/{it} -log {od}/efastats.log ; " \
                      "{me} -disperse {od}{it} -log {od}/disperse.log "\
        .format(me="muscle", od=outdir, it="/diversified_aln.efa", ou="diversified_aln_confseq.efa",
                ot="letterconf.afa", al="aln.fa", hl="letterconf.html", jv="letterconf_jalview.features")
    st_muscle_data = "Starting collection of additional muscle alignment data"
    ed_muscle_data = "Successfully obtained additional data"

    # Running analysis
    process_run(cmd_genbank, st_genbank, ed_genbank)
    process_run(cmd_blast_db, st_blast_db, ed_blast_db)
    process_run(cmd_blast_run, st_blast_run, ed_blast_run)
    multi_sanger()
    # process_run(cmd_muscle, st_muscle, ed_muscle)
    # process_run(cmd_muscle_data, st_muscle_data, ed_muscle_data)

    # Finishing analysis
    tend = int(time.time() - start_time)
    elapsed_time = "{:02d}:{:02d}:{:02d}".format(tend // 3600, (tend % 3600 // 60), tend % 60)
    logging.info("Analysis successful. Time to completion: {et}".format(et=elapsed_time))
    sys.exit(0)


if __name__ == '__main__':
    start_time = time.time()
    date = datetime.now()
    main()
