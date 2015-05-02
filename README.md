# banzai! #

This project analyzes the sequencing results from an Illumina MiSeq, run on a set of PCR amplicons.
It was written for a set of PCR amplicons which have been
The primary file is a Bash Shell script, which should run on Unix and Linux machines. The script makes heavy usage of Unix command line utilities (such as find, grep, sed, awk, and more) and is written for the BSD versions of those programs as found on standard installations of Mac OSX.

## Basic implementation ##
Simply edit the parameters in the file 'banzai_params.sh', then type into a terminal:
```sh
bash /Users/user_name/path/to/the/file/banzai.sh
```


## Dependencies ##
* seqtk # reverse complementing entire fastq/a files
* cutadapt # note that version 1.7 will support anchored 3' sequences.
* PEAR # merging paired-end reads
* usearch # quality filtering of merged paired-end fastq files, OTU clustering
* blast+ # taxonomic assignment
* MEGAN # taxonomic assignment
* R # ecological analyses. Requires the packages
* fastqc # quality control of raw sequencing fastq files

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
AATAGCGCTACGT ; READ1; READ2; READ3

# Note that the original script also ouput a file of the sequences only (no names), but I removed this functionality on 20150417
