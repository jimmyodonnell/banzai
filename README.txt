This project analyzes the sequencing results from an Illumina MiSeq, run on a sample of PCR amplicons.

The primary file is a Bash Shell script, which should run on Unix and Linux machines. The script makes heavy usage of Unix command line utilities (such as find, grep, sed, awk, and more) and is written for the BSD versions of those programs as found on standard installations of Mac OSX.

Dependencies:
seqtk
cutadapt # note that version 1.7 will support anchored 3' sequences.
PEAR
usearch
blast+
MEGAN
R

Wishlist/notes to self:
-(consider piping pear (and other/all?) output by adding > pear_log.txt to end of line)

- For now, this script will only work if the output of merged files in PEAR is less than 8GB. A little tweaking could fix this.

- Potential files to be removed as part of cleanup at the end of script:
homopolymer_line_numbers.txt
2_filtered.fasta
1_merged.assembled.fastq
1_merged.assembled_A.fastq
1_merged.assembled_B.fastq
