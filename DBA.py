#!/usr/bin/env python

import argparse
import logging
import multiprocessing
import os
import shutil
import subprocess
import sys
import time
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from datetime import datetime

# Credentials
__author__ = "M.D.C. Jansen"
__version__ = "1.0.2"
__date__ = "21/06/2023"


# Check thread input
def cpu_threads(thread_input):
    if multiprocessing.cpu_count() > thread_input:
        return thread_input
    else:
        return multiprocessing.cpu_count()


# Graceful preemptive termination of analysis
def terminate():
    tend = int(time.time() - start_time)
    elapsed_time = "{:02d}:{:02d}:{:02d}".format(tend // 3600, (tend % 3600 // 60), tend % 60)
    logging.error("Analysis terminated after: {et}\n\n\n".format(et=elapsed_time))
    sys.exit(1)


# Running analysis
def process_run(cmd_in, process_start, process_complete):
    logging.info(process_start)
    process = subprocess.run([cmd_in],
                             shell=True,
                             text=True,
                             capture_output=True
                             )
    if process.stderr is not None and len(process.stderr) != 0:
        if cmd_in == cmd_muscle or cmd_muscle_data:
            pass
        else:
            logging.error(process.stderr)
            terminate()
    else:
        logging.info(process_complete)


# Multithreading chromatogram production
def multi_sanger():
    processes = []
    for i in range(0, len(folder_names)):
        cmd = "sangerseq_viewer -s {wd}/{gbn}.gb -q {wd}/{fn} -o {od}/{cg}_{fn}.pdf -l 800" \
            .format(wd=workdir, gbn=genbank_nucl, fn=folder_names[int(i)], od=outdir, cg="chromatogram")
        st_log = "Creating chromatogram for {fn}".format(fn=folder_names[int(i)])
        ed_log = "Successfully produced chromatogram for {fn}".format(fn=folder_names[int(i)])
        p = multiprocessing.Process(target=process_run, args=(cmd, st_log, ed_log))
        processes.append(p)
        p.start()
    for j in processes:
        j.join()


# Reverse compliment input fasta
def rev_comp(record):
    return SeqRecord(seq=record.seq.reverse_complement(),
                     id="RC_{rid}".format(rid=record.id),
                     description="reverse complement"
                     )


def main():

    # Logging inputs
    logging.info("Settings barcoding analysis:\n\n"
                 "Input folder:\t\t\t\t{inp}\n"
                 "Input fasta:\t\t\t\t{inf}\n"
                 "Input mao:\t\t\t\t{inm}\n"
                 "Output directory:\t\t\t{od}\n"
                 "Genbank nucleotide accession no.:\t{gbn},\t{gba}\n"
                 "Genbank protein accession no.:\t\t{gbp}\n"
                 "Threads:\t\t\t\t{td}\n"
                 "Reverse complement fasta:\t\t{rev}\n"
                 "Keep all files:\t\t\t\t{kp}\n"
                 .format(inp=argument.i, inf=input_fasta, inm=input_mao, od=os.path.basename(outdir), gba=genbank_out,
                         gbn=genbank_nucl, gbp=genbank_prot, td=threads, rev=rev, kp=keep))

    # Analysis functions
    # Obtaining data from genbank
    cmd_genbank = "bio fetch {gbn} > {wd}/{gbn}.gb ; " \
                  "bio fetch {gbn} --format fasta > {wd}/{gbn}.fasta ; " \
                  "bio fetch {gbp} --format fasta > {wd}/{gbp}.prot ; " \
                  "bio fetch {gba} --format fasta > {wd}/{gba}.fasta"\
        .format(wd=workdir, gbn=genbank_nucl, gbp=genbank_prot, gba=genbank_out)
    st_genbank = "Obtaining genbank entries:\t\t\t{gbn},\t{gbp},\t{gba}"\
        .format(gbn=genbank_nucl, gbp=genbank_prot, gba=genbank_out)
    ed_genbank = "Successfully obtained Genbank entries:\t{gbn},\t{gbp},\t{gba}"\
        .format(gbn=genbank_nucl, gbp=genbank_prot, gba=genbank_out)

    # Running BLASTn and BLASTx
    cmd_blast_db = "makeblastdb -in {wd}/{gbn}.fasta -dbtype nucl -out {wd}/{gbn}_nucldb ; " \
                   "makeblastdb -in {wd}/{gbp}.prot -dbtype prot -out {wd}/{gbp}_protdb"\
        .format(wd=workdir, gbn=genbank_nucl, gbp=genbank_prot)
    st_blast_db = "Creating blast databases"
    ed_blast_db = "Successfully created databases"
    cmd_blast_run = "blastn -query {wd}/{inf} -db {wd}/{gbn}_nucldb -out {od}/blastn.tsv -outfmt 6 ; " \
                    "blastx -query {wd}/{inf} -db {wd}/{gbp}_protdb -out {od}/blastx.tsv -outfmt 6"\
        .format(wd=workdir, inf=input_fasta, od=outdir, gbn=genbank_nucl, gbp=genbank_prot)
    st_blast_run = "Starting BLASTn and BLASTx analysis"
    ed_blast_run = "Successfully completed BLASTn and BLASTx analysis"

    # Aligning sequences with MUSCLE
    global cmd_muscle, cmd_muscle_data
    cmd_muscle = "muscle -align {wd}/{sq} -output {od}/{ou} -threads {ts} && " \
                 "muscle -align {wd}/{sq} -output {od}/{ot} -threads {ts} -diversified" \
        .format(wd=workdir, sq="aligning.fasta", od=outdir, ou="aln.fa", ot="diversified_aln.fa", ts=threads)
    st_muscle = "Starting alignment with muscle"
    ed_muscle = "Successfully performed alignment"
    cmd_muscle_data = "{me} -addconfseq {od}/{it} -output {od}/{odc} ; " \
                      "{me} -letterconf {od}/{it} -ref {od}/{al} -html {od}/{hl}" \
        .format(me="muscle", od=outdir, al="aln.fa", oac="aln_confseq.fa", it="diversified_aln.fa",
                ot="letterconf.nfa", odc="diversified_confseq.fa", hl="letterconf.html")
    st_muscle_data = "Starting collection of additional muscle alignment data"
    ed_muscle_data = "Successfully obtained additional data"
    cmd_muscle_max = "muscle -maxcc {od}/{ot} -output {od}/maxcc_{ot}".format(od=outdir, ot="diversified_aln.fa")
    st_muscle_max = "Starting extraction of maximum confidence alignment"
    ed_muscle_max = "Successfully extracted alignment with highest confidence"

    # Utilising Jalview for manual review and trimming of alignment
    cmd_jalview = "jalview -open {od}/{imc}".format(od=outdir, imc="maxcc_diversified_aln.fa")
    st_jalview = "Opening Jalview for manual trimming"
    ed_jalview = "Jalview closed, continuing analysis"

    # Constructing neighbour-joining tree with MEGA
    cmd_mega = "megacc -a {wd}/{mao} -d {od}/{efa} -o {od}"\
        .format(wd=workdir, mao=input_mao, od=outdir, efa="maxcc_diversified_aln.fa")
    st_mega = "Starting phylogenetic analysis with MEGA X"
    ed_mega = "Successfully performed phylogenetic analysis"

    # Running analysis
    process_run(cmd_genbank, st_genbank, ed_genbank)
    multi_sanger()
    process_run(cmd_blast_db, st_blast_db, ed_blast_db)
    process_run(cmd_blast_run, st_blast_run, ed_blast_run)
    os.system("cat {wd}/{gbn}.fasta {wd}/{gba}.fasta {wd}/{inf} > {wd}/aligning.fasta"
              .format(wd=workdir, inf=input_fasta, gbn=genbank_nucl, gba=genbank_out))
    process_run(cmd_muscle, st_muscle, ed_muscle)
    process_run(cmd_muscle_data, st_muscle_data, ed_muscle_data)
    process_run(cmd_muscle_max, st_muscle_max, ed_muscle_max)
    process_run(cmd_jalview, st_jalview, ed_jalview)
    process_run(cmd_mega, st_mega, ed_mega)

    # Finishing analysis
    if keep is True:
        logging.info("Moving working directory to output")
        shutil.move(workdir, outdir)
    else:
        logging.info("Clearing working directory")
        shutil.rmtree(workdir)
    tend = int(time.time() - start_time)
    elapsed_time = "{:02d}:{:02d}:{:02d}".format(tend // 3600, (tend % 3600 // 60), tend % 60)
    logging.info("Successfully completed barcoding analysis\n\t\t\t\t     "
                 "Analysis completed after: {et}\n\n\n".format(et=elapsed_time))
    sys.exit(0)


if __name__ == '__main__':
    start_time = time.time()
    date = datetime.now()

    # Setup parser
    parser = argparse.ArgumentParser(prog="DBA",
                                     description="DBA is designed to automate the process of DNA barcoding "
                                                 "by utilising standard Sanger sequencing data "
                                                 "and user provided reference gene and protein accession numbers.\n"
                                                 "This pipeline obtains nucleic and protein sequences from genbank, "
                                                 "followed by chromatogram production in .pdf format by utilising "
                                                 "sangerseq viewer.\nNext, BLASTx and BLASTn are run to assess query "
                                                 "coverage against the reference gene.\nMUSCLE is used to align the "
                                                 "sequences, followed by manual review in jalview and "
                                                 "phylogenetic analysis by MEGA.",
                                     usage="%(prog)s -i <inputfolder> -n <genbank reference NC_ID> "
                                           "-y <genbank reference YP_ID> -g <genbank outgroup NC_ID> [options]",
                                     epilog="",
                                     formatter_class=argparse.RawTextHelpFormatter
                                     )
    parser._optionals.title = "List of arguments"
    parser.add_argument("-v", "--version",
                        help="Prints program version",
                        action="version",
                        version="Version: {v} Date: {d} By {a}".format(v=__version__, d=__date__, a=__author__))
    parser.add_argument("-i",
                        metavar="[input folder]",
                        help="Input folder containing per species folders, "
                             "a single concatenated fasta file for analysis and "
                             "a MEGA mao file for phylogenetic analysis",
                        required=len(sys.argv) != 1)
    parser.add_argument("-n",
                        metavar="[genbank reference nucl acc no.]",
                        help="genbank nucleotide accession number for the gene to be used as reference, "
                             "usually starts with the identifier NC_",
                        required=len(sys.argv) != 1)
    parser.add_argument("-y",
                        metavar="[genbank reference prot acc no.]",
                        help="genbank protein accession number for the protein to be used as reference,"
                             " usually starts with the identifier YC_",
                        required=len(sys.argv) != 1)
    parser.add_argument("-g",
                        metavar="[genbank nucl acc no. for out group]",
                        help="genbank nucleotide accession number for outgroup gene used during phylogenetic analysis, "
                             "usually starts with the identifier NC_",
                        required=len(sys.argv) != 1)
    parser.add_argument("-o",
                        metavar="[output folder]",
                        help="Output directory. "
                             "\nA default output folder will be produced in the following format will be created "
                             "if it hasn't been specified:"
                             "\nbarcoding_output_current-date_current-time",
                        required=False,
                        default="barcoding_output_{dt}/".format(dt=date.strftime("%d-%m-%Y_%H-%M-%S")))
    parser.add_argument("-t",
                        metavar="[cpu threads]",
                        help="Maximum amount of threads to be utilised during analysis."
                             "\nDefault: 20",
                        required=False,
                        default=cpu_threads(20),
                        type=int)
    parser.add_argument("-keep",
                        metavar='',
                        help="Keep all files produced during analysis."
                             "\nFiles are stored in 'workdir' folder at the location where this script"
                             "has been executed."
                             "\nDefault: False",
                        required=False,
                        default=False,
                        type=bool,
                        nargs="?",
                        const=True)
    parser.add_argument("-rev",
                        metavar='',
                        help="Reverse compliment the fasta input."
                             "\nReverse complimented file is saved as an additional file."
                             "\nDefault: False",
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
    genbank_nucl = str(argument.n)
    genbank_prot = str(argument.y)
    genbank_out = str(argument.g)
    threads = str(argument.t)
    keep = argument.keep
    rev = argument.rev

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

    # Validating inputfolder
    inputfolder = os.path.abspath(argument.i)
    if len(inputfolder) != 0 and os.path.exists(inputfolder):
        logging.info("Found folder:\t\t\t\t{fn}".format(fn=inputfolder))
        folder_names = [dir for dir in os.listdir(inputfolder) if os.path.isdir(os.path.join(inputfolder, dir))]
        folder_count = len(folder_names)
    else:
        logging.error("No input folder found.")
        terminate()

    # Defining input files
    suffix = [".fasta", ".mao"]
    for s in range(0, len(suffix)):
        input_file = [file for file in os.listdir(inputfolder) if file.endswith("{end}".format(end=suffix[int(s)]))
                      if os.path.isfile(os.path.join(inputfolder, file))]
        if len(input_file) == 1:
            logging.info("Found {end} file:\t\t\t\t{inf}"
                         .format(end=suffix[int(s)], inf=os.path.join(inputfolder, input_file[0])))
            if input_file[0].endswith(".fasta"):
                input_fasta = input_file[0]
            elif input_file[0].endswith(".mao"):
                input_mao = input_file[0]
        else:
            logging.error("Incorrect number of input {end} files found\n\t\t\t\t     "
                          "Expected one, found {nf}".format(end=suffix[int(s)], nf=len(input_file)))
            terminate()

    # Validating output
    outdir = os.path.abspath(argument.o)
    if os.path.exists(outdir):
        logging.warning("Output directory already exists\n\t\t\t\t     Continuing analysis in default output folder")
        outdir = "barcoding_output_{dt}/".format(dt=date.strftime("%d-%m-%Y_%H-%M-%S"))
        logging.info("Creating output directory:\t\t\t{od}".format(od=outdir))
        os.makedirs(outdir)
    else:
        logging.info("Creating output directory:\t\t\t{od}".format(od=outdir))
        os.makedirs(outdir)

    # Preparing working environment
    logging.info("Preparing working folder for {fc} inputs"
                 .format(fc=str(folder_count)))
    workdir = os.path.join(root, "workdir")
    if not os.path.isdir(workdir):
        os.makedirs(workdir)
    else:
        shutil.rmtree(workdir)
        os.makedirs(workdir)
    num_folder = 0
    while num_folder != folder_count:
        os.makedirs("{wd}/{fn}"
                    .format(wd=workdir, fn=folder_names[num_folder]))
        os.system("cp {inp}/{fn}/*/* {wd}/{fn}/"
                  .format(inp=inputfolder, fn=folder_names[num_folder], wd=workdir))
        num_folder += 1
    shutil.copy(os.path.join(inputfolder, input_fasta), workdir)
    shutil.copy(os.path.join(inputfolder, input_mao), workdir)
    shutil.move(logfile, outdir)
    logfile = "{od}/barcoding.log".format(od=outdir)
    if rev is True:
        logging.info("Starting reverse complementing:\t\t{nf}".format(nf=input_fasta))
        rev_fasta = map(rev_comp, SeqIO.parse("{wd}/{inf}".format(wd=workdir, inf=input_fasta), "fasta"))
        SeqIO.write(rev_fasta, "{wd}/rev_comp_{inf}".format(wd=workdir, inf=input_fasta), "fasta")
        input_fasta = "rev_comp_{inf}".format(inf=input_fasta)
        logging.info("Successfully reverse complemented to:\t{nf}".format(nf=input_fasta))

    # Start analysis
    main()
