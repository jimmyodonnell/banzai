#!/bin/bash

# What is the path to the primer tags?
# This file should be simply a list of sequences, one per line, of each of the tags, WITH A TRAILING NEWLINE!
# To make a trailing newline, make sure when you open the file, you have hit enter after the final sequence.
PRIMER_TAGS='/Users/threeprime/Documents/Data/IlluminaData/12S/12S_Tags.txt'
# What is the maximum number of Ns to allow at the end of a sequence before a tag is reached? NOTE: currently, this will only matter for the tag on the 3' end.
# TAG_N_MAX="9" # THIS IS NOT WORKING YET. SET TO DEFAULT 9

# What is the path to the reads?
READ1='/Users/threeprime/Documents/Data/IlluminaData/12S/20140930/PCpoolC_S1_L001_R1_001.fastq.gz'
READ2='/Users/threeprime/Documents/Data/IlluminaData/12S/20140930/PCpoolC_S1_L001_R2_001.fastq.gz'

# Is it ok to rename the sequences within a fasta file
# This will happen after the fastq has been converted to a fasta file at the quality filtering step.
RENAME_READS="YES"

# What is the maximum expected length of the fragment of interest, including primers?
LENGTH_FRAG="160"

# What is the length of the reads of the Illumina run? (i.e. how long are the sequences in each of the run fastq files (R1 and R2)?)
LENGTH_READ="150"

# Specify the path to the MEGAN executable file you want to use.
megan_exec='/Applications/megan/MEGAN'

# PRIMER REMOVAL
# Specify a path to the fasta file containing the two primers used to generate the amplicons you sequenced:
PRIMER_FILE='/Users/threeprime/Documents/Data/IlluminaData/12S/primers_12S.fasta'
PRIMER_MISMATCH_PROPORTION="0.10"

# Would you like to remove reads containing runs of consecutive identical bases (homopolymers)?
REMOVE_HOMOPOLYMERS="NO"
# What is the maximum homopolymer length you're willing to accept?
# Reads containing runs of identical bases longer than this will be discarded.
HOMOPOLYMER_MAX="7"

# CLUSTERING:
# Note that this must be an integer. Contact me if this is a problem
CLUSTERING_PERCENT="99"

# If you'd like to blast reads without clustering them first, set this to YES:
BLAST_WITHOUT_CLUSTERING="YES"

# BLAST:
# Specify the path to the BLAST database.
# Note this should be a path to any one of three files WITHOUT their extension *.nhr, *.nin, or *.nsq
BLAST_DB='/Users/threeprime/Documents/Data/genbank/12S/Vert12sdb/Vert12S_BLAST_DB/Vert_full_mtDNA_genomes_RefSeq+Bony_Cart_Cetacean_12S_partial'
# BLAST PARAMETERS
PERCENT_IDENTITY="95"
WORD_SIZE="80"
EVALUE="0.001"
# number of matches recorded in the alignment:
MAXIMUM_MATCHES="5"

# What is the lowest taxonomic rank at which MEGAN should group OTUs?
COLLAPSE_RANK="Genus"
MINIMUM_SUPPORT="2"
MINIMUM_COMPLEXITY="0"
TOP_PERCENT="2"
MINIMUM_SUPPORT_PERCENT="0"
MINIMUM_SCORE="150"

# Would you like to delete extraneous intermediate files once the analysis is finished? YES/NO
PERFORM_CLEANUP="YES"


####################### WOULD YOU LIKE TO PICK UP FROM AN EXISTING FILE?
# If reanalyzing existing demultiplexed data, point this variable to the directory storing the individual tag folders.
EXISTING_DEMULTIPLEXED_DIR='/Users/threeprime/Documents/Data/IlluminaData/16S/20141020/Analysis_20141023_1328/demultiplexed'

# Have the reads already been paired?
ALREADY_PEARED="YES"
# YES/NO
PEAR_OUTPUT='/Users/threeprime/Documents/Data/IlluminaData/12S/20140930/Analysis_20141030_2020/1_merged.assembled.fastq'

# Have the merged reads been quality filtered?
ALREADY_FILTERED="YES" # YES/NO
FILTERED_OUTPUT='/Users/threeprime/Documents/Data/IlluminaData/12S/20140930/Analysis_20141030_2020/2_filtered.fasta'

# Should demultiplexed samples be concatenated for annotation as a single unit? (Each read can still be mapped back to samples)
CONCATENATE_SAMPLES="YES"
