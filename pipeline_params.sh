#!/bin/bash

# An attempt to cause the script to exit if any of the commands returns a non-zero status (i.e. FAILS).
set -e

# Are the reads multiplexed with tagged primers? (YES/NO)
MULTIPLEXED="YES"

# What is the path to the primer tags?
# This file should be simply a list of sequences, one per line, of each of the tags, WITH A TRAILING NEWLINE!
# To make a trailing newline, make sure when you open the file, you have hit enter after the final sequence.
PRIMER_TAGS='/Users/threeprime/Documents/Data/IlluminaData/12S/12s_Tags.txt'

# IF DATA NEEDS TO BE DEMULTIPLEXED FIRST:
# What is the path to the reads?
READ1='/Users/threeprime/Documents/Data/IlluminaData/12S/20140930/PCpoolC_S1_L001_R1_001.fastq'
READ2='/Users/threeprime/Documents/Data/IlluminaData/12S/20140930/PCpoolC_S1_L001_R2_001.fastq'

# What is the maximum expected length of the fragment of interest, including primers?
LENGTH_FRAG="162"

# What is the length of the reads of the Illumina run? (i.e. how long are the sequences in each of the run fastq files (R1 and R2)?)
LENGTH_READ="150"

# Specify the path to the MEGAN executable file you want to use.
megan_exec='/Applications/megan/MEGAN'

# PRIMER REMOVAL
# Specify a path to the fasta file containing the two primers used to generate the amplicons you sequenced:
PRIMER_FILE='/Users/threeprime/Documents/Data/IlluminaData/16S/primers_16S.fasta'
PRIMER_MISMATCH_PROPORTION="0.10"

# What is the maximum number of consecutive identical bases you're willing to accept?
# Reads containing runs of identical bases longer than this will be discarded.
HOMOPOLYMER_MAX="7"

# CLUSTERING:
# Note that this must be an integer. Contact me if this is a problem
CLUSTERING_PERCENT="97"

# BLAST:
# Specify the path to the BLAST database.
# Note this should be a path to any one of three files WITHOUT their extension *.nhr, *.nin, or *.nsq
BLAST_DB='/Users/threeprime/Documents/Data/genbank/16Smetazoa/16Smetazoa'
# BLAST PARAMETERS
PERCENT_IDENTITY="95"
WORD_SIZE="100"
EVALUE="10"
# number of matches recorded in the alignment:
MAXIMUM_MATCHES="20"

# At what taxonomic rank should MEGAN group OTUs
COLLAPSE_RANK="Family"

# Would you like to delete extraneous intermediate files once the analysis is finished? YES/NO
PERFORM_CLEANUP="NO"
