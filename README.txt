This project analyzes the sequencing results from an Illumina MiSeq, run on a sample of PCR amplicons.

The primary file is a Bash Shell script, which should run on Unix and Linux machines. The script makes heavy usage of Unix command line utilities (such as find, grep, sed, awk, and more) and is written for the BSD versions of those programs as found on standard installations of Mac OSX.

Dependencies:
seqtk # reverse complementing entire fastq/a files; perhaps unnecessary now.
cutadapt # note that version 1.7 will support anchored 3' sequences.
PEAR # merging paired-end reads
usearch # various tasks, including OTU clustering
blast+ # taxonomic assignment
MEGAN # taxonomic assignment
R # ecological analyses
fastqc # quality control of raw sequencing fastq files

Wishlist/notes to self:
-(consider piping pear (and other/all?) output by adding > pear_log.txt to end of line)

- For now, this script will only work if the output of merged files in PEAR is less than 8GB. A little tweaking could fix this.

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
