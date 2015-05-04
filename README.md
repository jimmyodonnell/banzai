#banzai!#
=========

**banzai** is a BASH shell script that links together the disparate analyses needed to process the raw sequencing results from an Illumina run, the eventual goal being a table of the number of sequences per taxon found in a set of samples. Preliminary ecological analyses are included as well.

The script should run on Unix and Linux machines. The script makes heavy usage of Unix command line utilities (such as find, grep, sed, awk, and more) and is written for the BSD versions of those programs as found on standard installations of Mac OSX. I tried to use POSIX-compliant commands wherever possible.

## Basic implementation ##
Simply edit the parameters in the file 'banzai_params.sh', then type into a terminal:

```sh
bash /Users/user_name/path/to/the/file/banzai.sh
```


## Dependencies ##
Aside from the standard command line utilities that are already included on Unix machines (awk, sed, grep, etc), this script relies on the following tools:

* **PEAR**: merging paired-end reads
* **usearch**: quality filtering of merged paired-end fastq files, OTU clustering
* **cutadapt**: primer removal
* **seqtk**: reverse complementing entire fastq/a files
* **python**: fast consolidation of duplicate sequences (installed by default on Macs)
* **blast+**: taxonomic assignment
* **MEGAN**: taxonomic assignment
* **R**: ecological analyses. Requires the packages **vegan** and **gtools**
* **fastqc**: quality control of raw sequencing fastq files

Wishlist/TODO/notes to self:
- streamline config file

- Potential files to be removed as part of cleanup at the end of script:
homopolymer_line_numbers.txt
2_filtered.fasta
1_merged.assembled.fastq
1_merged.assembled_A.fastq
1_merged.assembled_B.fastq

- 2014 11 12: Noticed that the reverse tag removal step removes the tag label from the sequenceID line of fasta files if the tag sequence is RC-palindromic!

ADD LIBRARY NAMES VARIABLE

###############
# LOG FILE
###############
# logging implemented! But...

INVOKE THE SCRIPT BY TYPING 'bash script.sh' NOT 'sh script.sh' !!!!!


# Alternatively, to have the pipeline print to terminal AND file:
sh script.sh  2>&1 | tee ~/Desktop/logfile.txt



# Removal of duplicate sequences (dereplicate_fasta.py)
Input: a fasta file (e.g. 'infile.fasta')
Output: a file with the same name as the input but with the added extension '.all' (e.g. 'infile.fasta.all')
This file contains each unique DNA sequence from the fasta file, followed by the labels of the reads matching this sequence
Thus, if an input fasta file consisted of three reads with identical DNA sequences:

	>READ1
	AATAGCGCTACGT
	>READ2
	AATAGCGCTACGT
	>READ3
	AATAGCGCTACGT

The output file is as follows:

	AATAGCGCTACGT; READ1; READ2; READ3

# Note that the original script also ouput a file of the sequences only (no names), but I removed this functionality on 20150417
