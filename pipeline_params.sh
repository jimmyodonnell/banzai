#!/bin/bash

# YOUR analysis folder must contain:
# analysis_pipeline.sh
# pear_params.sh
# primers.fasta # A fasta file containing only the forward and reverse primers, in that order, in their 5'-3' orientation.


# What is the maximum expected length of the fragment of interest, including primers?
LENGTH_FRAG="162"

# What is the length of the reads of the Illumina run? (i.e. how long are the sequences in each of the run fastq files (R1 and R2)?)
LENGTH_READ="150"

# Specify the path to the MEGAN executable file you want to use.
megan_exec='/Applications/megan/MEGAN'

# Specify the path to the BLAST database.
# Note this should be a path to any one of three files WITHOUT their extension *.nhr, *.nin, or *.nsq
BLAST_DB='/Users/threeprime/Documents/Data/genbank/16Smetazoa/16Smetazoa'

# PRIMER REMOVAL
# Specify a path to the fasta file containing the two primers used to generate the amplicons you sequenced:
PRIMER_FILE='/Users/threeprime/Documents/Data/IlluminaData/16S/primers_16S.fasta'
PRIMER_MISMATCH_PROPORTION="0.10"

# What is the maximum number of consecutive identical bases you're willing to accept?
# Reads containing runs of identical bases longer than this will be discarded.
HOMOPOLYMER_MAX="7"

#CLUSTERING:
CLUSTERING_PERCENT="97"

#BLAST:
PERCENT_IDENTITY="95"
WORD_SIZE="100"
# number of matches recorded in the alignment:
MAXIMUM_MATCHES="20"

# Would you like to delete extraneous intermediate files once the analysis is finished? YES/NO
PERFORM_CLEANUP="YES"
