# DBA

DNA barcoding analysis pipeline designed to perform phylogenetic analysis on raw sanger sequencing data.

## Table of contents

* [Requirements](#Requirements-and-dependencies)
* [Installation](#Installation)
* [Example commands](#Example-commands)
* [Full usage](#Full-usage)
* [Method](#Method)
* [Interpeting output](#Interpeting-output)
* [Acknowledgements](#Acknowledgements)
* [License](#License)


## Requirements and dependencies
DBA is designed to run on linux-based systems. The pipeline was specifcally tested on Ubuntu 22.04.2 LTS.
Some tools and most dependencies used during analysis will be installed automatically when the conda environment is being constructed.

The following prerequisites should be installed prior to using DBA
* Anaconda3
* Jalview
	- Note: Jalview should be a system-wide executable as "jalview"
* Java	
* MEGA11
* Python 3.8+
* Python 3.8+
* Java
* Jalview
	- Note: Jalview should be a system-wide executable as "jalview"
* MEGA11
* Anaconda3


## Installation

#### Install dependecies:
```

```

#### Install DBA:
```
git clone https://github.com/mdcjansen/DBA
cd path/to/DBA
conda env create -f environment.yml
chmod a+x DBA.py
```

#### Allow DBA to be executed from anywhere within the anaconda environment:
```
ln -s /path/to/DBA.py /path/to/anaconda3/envs/DBA/bin/DBA
```


## Example commands

#### Activate conda environment for analysis:
```
conda activate DBA
```

#### Display options:
This command will display all required and optional arguments in the format presented [below](#Full-usage)
```
DBA.py
DBA.py -h "or" --help
```

#### Minimum required input:
```
DBA.py -i <input folder> -n <genbank reference NC_ID> -y <genbank reference YC_ID> -g <genbank outgroup NC_ID>
```

#### Example input:
The input below was used to produce the case data available within this repository:
```
DBA.py -i case_data -n NC_009065 -y YP_001054869 -g NC_02842
```

## Full usage
```
usage: DBA -i <inputfolder> -n <genbank reference NC_ID> -y <genbank reference YP_ID> -g <genbank outgroup NC_ID> [options]

DBA is designed to automate the process of DNA barcoding by utilising standard Sanger sequencing data and user provided reference gene and protein accession numbers.
This pipeline obtains nucleic and protein sequences from genbank, followed by chromatogram production in .pdf format by utilising sangerseq viewer.
Next, BLASTx and BLASTn are run to assess query coverage against the reference gene.
MUSCLE is used to align the sequences, followed by manual review in jalview and phylogenetic analysis by MEGA.

List of arguments:
  -h, --help            show this help message and exit
  -v, --version         Prints program version
  -i [input folder]     Input folder containing per species folders, a single concatenated fasta file for analysis and a MEGA mao file for phylogenetic analysis
  -n [genbank reference nucl acc no.]
                        genbank nucleotide accession number for the gene to be used as reference, usually starts with the identifier NC_
  -y [genbank reference prot acc no.]
                        genbank protein accession number for the protein to be used as reference, usually starts with the identifier YC_
  -g [genbank nucl acc no. for out group]
                        genbank nucleotide accession number for outgroup gene used during phylogenetic analysis, usually starts with the identifier NC_
  -o [output folder]    Output directory. 
                        A default output folder will be produced in the following format will be created if it hasn't been specified:
                        barcoding_output_current-date_current-time
  -t [cpu threads]      Maximum amount of threads to be utilised during analysis.
                        Default: 20
  -keep []              Keep all files produced during analysis.
                        Files are stored in 'workdir' folder at the location where this scripthas been executed.
                        Default: False
  -rev []               Reverse compliment the fasta input.
                        Reverse complimented file is saved as an additional file.
                        Default: False
```


## Method
DBA requires a specific folder structure in order to perform the analysis.
A visualisation of this folder structure is given below where case_data is given as `-i` argument:

```
	.
	├── ...
	├── case_data               
	│	├── sanger.fasta				# FASTA file containing all sanger sequences to be used during analysis
	│	├── phylo.mao					# MEGA11 generated file used to perform phylogenetic analysis
	│	├── species_01					# Folder with sequencing data of first species to be analysed
	│	│	├── sample_01				# Folder containing sanger sequencing data of the first sample
	│	│	│	├── sanger_data.ab1
	│	│	│	└── ...
	│	│	└── ...
	│	├── species_02					# Folder with sequencing data of second species to be analysed
	│	│	├── sample_01				# Folder containing sanger sequencing data of the first sample
	│	│	│	├── sanger_data.ab1
	│	│	│	└── ...
	│	│	└── ...
	│	└── ...
	└── ...
```

DBA will perform phylogenetic analysis as presented below:
1. If specified, DBA will reverse complement the input fasta file
2. DBA will obtain genbank information on input genbank accession numbers
3. [sangerseq_viewer](https://github.com/ponnhide/sangerseq_viewer) is utilised to produce a chromatogram for each of the input species found
4. [BLAST](https://doi.org/10.1016/S0022-2836(05)80360-2) is performed on the input fasta sequence against the input reference sequences
5. Multple sequence alignment is performed with [MUSCLE](https://doi.org/10.1101/2021.06.20.449169)
6. DBA initates [Jalview](https://doi.org10.1093/bioinformatics/btp033) for manual review and trimming of alignment data
7. [MEGA11](https://doi.org/10.1093/molbev/msab120) is used to perform phylogentic analysis and produce a consensus tree
8. If specified, supplemental data used during analysis is moved to the output folder and the analysis ends.


## Interpeting output
DBA produces a single output folder where all results can be found.
All procedures during analysis are logged and stored in the barcoding.log file

#### BLAST output
BLAST results are placed into two tsv files named `blastn.tsv` and `blastx.tsv` for BLASTn and BLASTx results respectively.
Standard BLAST tabular output format 6 is used to generate the output.
An example with case data and added header is presented below:

qseqid|sseqid|pident|length|mismatch|gapopen|qstart|qend|sstart|send|evalue|bitscore
---|---|---|---|---|---|---|---|---|---|---|---|
RC_EF71414350_EF71414350|NC_009065.1|98.214|672|8|4|13|680|5513|6184|0.0|1171
RC_EF71414351_EF71414351|NC_009065.1|99.519|624|3|0|5|628|5507|6130|0.0|1136
RC_EF71414352_EF71414352|NC_009065.1|97.765|671|9|6|5|670|5507|6176|0.0|1151
RC_EF71414353_EF71414353|NC_009065.1|98.336|661|6|5|10|670|5513|6168|0.0|1155
RC_EF71414354_EF71414354|NC_009065.1|99.522|628|3|0|6|633|5507|6134|0.0|1144

#### Sangerseq_viewer output
Sangerseq_viewer produces a chromatogram in pdf format for each species folder present in the input folder.
These chromatograms display all sequences found inside a particular species folder and present mapped genbank genes.
An example produced from case data is presented below:

![Screenshot](chromatogram.png)

#### MUSCLE output
MUSCLE produces multiple alignments files, all presented with the .fa suffix.
The alignment with the highest column confidence is extracted from the diversified ensemble and opened in jalview to be reviewed by the user.
Additionally, an html file will be produced for the diversified ensemble containing the letter confidence values of the ensemble.

#### MEGA output
MEGA generates phylogenetic trees in newick format using user specified parameters specified within the mao file.
Newick trees and consensus tree are both present in the output folder. Summary and partition text files are also saved within the output.


## Acknowledgements
DBA uses the following tools in the pipeline:
* [BLAST](https://doi.org/10.1016/S0022-2836(05)80360-2)
* [Jalview](https://doi.org10.1093/bioinformatics/btp033)
* [MEGA11](https://doi.org/10.1093/molbev/msab120)
* [MUSCLE](https://doi.org/10.1101/2021.06.20.449169)
* [sangerseq_viewer](https://github.com/ponnhide/sangerseq_viewer)


## License
[MIT license](https://opensource.org/license/mit/)

