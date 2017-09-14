#!/usr/bin/env bash


################################################################################
# INPUT
################################################################################
# What is the file path to the directory containing all of the libraries/reads?
PARENT_DIR="${BANZAI_DIR}"/data/PS_16S

# Where is the sequencing metadata file? (SEE FORMATTING GUIDELINES IN README!)
SEQUENCING_METADATA="${PARENT_DIR}"/metadata.csv


################################################################################
# OUTPUT
################################################################################
# This script will generate a directory (folder) containing the output of the script.
# Where do you want this new folder to go?
OUTPUT_DIRECTORY="${HOME}" #"${PARENT_DIR%/*}"

# You can optionally specify a folder into which the script copies a PDF containing some results.
# The pdf is created by default in the analysis folder specified above, but
# if you set this to your DropBox or Google Drive Folder, you can check it out from anywhere.
OUTPUT_PDF_DIR=""


################################################################################
# METADATA DETAILS
################################################################################
# Specify columns for raw sequencing files:
COLNAME_FILE1="file1"
COLNAME_FILE2="file2"

# MUST be unique for each row!
COLNAME_SAMPLE_ID="sample_id"

COLNAME_DNA_ID="dna_id"
COLNAME_SAMPLE_TYPE="dna_source"

# Your metadata must have a column corresponding to the subfolders containing the raw reads.
# In order to make this flexible across both multiple and single library preps, you must include this even if you only sequenced one library (sorry!).
COLNAME_ID1_NAME="pri_index_name"
COLNAME_ID1_SEQ="pri_index_seq"

################################################################################
# MERGE PAIRED READS
################################################################################
USE_PEAR_DEFAULTS="NO"
# For more information on these parameters, type into a terminal window: pear -help
# Bokulich recommends:
# Quality_Threshold=3, r=3 (PEAR only considers r=2), UNCALLEDMAX=0
# TRIMMIN= 0.75 * LENGTH_READ # this is hard-coded in the script banzai.sh

# TODO grab this from a fragment_size column in the sequencing metadata file
### ***** REMEMBER TO WATCH FOR ZEROS WHEN IMPLEMENTING THIS!
# Is there a column in the metadata file for fragment size?
insert_size_in_metadata="NO"
# If YES, what is the name?
COLNAME_INSERT_SIZE="insert_size"

# do you want banzai to automatically calculate the expected assembled sequence lengths and overlap based on read length and fragment size?
calculate_merge_length="NO" # [ YES | NO]

# If fragment size is not specified in metadata, specify it here.
# What is the maximum expected length of the fragment of interest?
# This is the length of the fragments input into library prep --
# i.e. with (indexed) primers, but without library index or sequencing adapters
LENGTH_FRAG="180"

# if "NO", provide the following values for PEAR:
minimum_overlap="10" # [10]
assembled_max="10000" # [1000]
assembled_min="50" # [50]
# note that as of 20151124, the PEAR default for assembly max (0) doesn't merge any reads.

# --quality-threshold
Quality_Threshold=15

# proportion of allowed uncalled bases (--max-uncalled-base)
UNCALLEDMAX=0

# which statistical test (--test-method)
TEST=1

# cutoff p-value (--p-value)
PVALUE=0.01

# scoring method type (--score-method)
SCORING=2

# What is the minimum final sequence length you'd like to include in analyses?
# Bokulich et al. recommend 75% of the expected fragment size.
# this includes PCR primers, but not anything ligated on during library prep.
min_seq_length=75
# equivalent variable was TRIMMIN

################################################################################
# QUALITY FILTERING
################################################################################
# Substantial quality filtering (e.g. trimming, minimum length, etc) is performed by PEAR during read merging.
# You may also want to exclude sequences containing more than a specified threshold of 'expected errors'
# This number is equal to the sum of the error probabilities.
# For more information on this parameter, Google the usearch help
Perform_Expected_Error_Filter="YES" # [YES|NO]
Max_Expected_Errors="0.5"

################################################################################
# HOMOPOLYMERS
################################################################################
# 454 sequencers have trouble correctly identifying the length of a strech of consecutive identical bases (homopolymers).
# Illumina machines do not have this problem, yet some people are still paranoid.
# Would you like to remove reads containing runs of consecutive identical bases (homopolymers)?
REMOVE_HOMOPOLYMERS="NO"
# What is the maximum homopolymer length you're willing to accept?
# Reads containing runs of identical bases longer than this will be discarded.
HOMOPOLYMER_MAX="10"


################################################################################
# DEMULTIPLEXING
################################################################################

# Do the reads contain index sequences which identifies their sample of origin?
SECONDARY_INDEX="YES"

# Specify the nucleotide sequences that differentiate multiplexed samples
# (sometimes, confusingly referred to as "tags" or "barcodes")
# these are the secondary index -- the primary index added with the sequencing adapters should not be in the sequence data
# You can grab these from the file specified above (SEQUENCING_METADATA) by specifying the column name of index sequences.
COLNAME_ID2_SEQ="sec_index_seq"

# How many nucleotides pad the 5' end of the tag sequence?
# TODO build in flexibility (this number is unused right now)
TAG_Ns="3"
SECONDARY_INDEX_START="4"
COLNAME_ID2_START="sec_index_start"

################################################################################
# PRIMER REMOVAL
################################################################################
# Specify the primers used to generate these amplicons.
# As with the multiplex indexes, Banzai will grab these from the file SEQUENCING_METADATA.
# You must indicate the column names of the forward and reverse primers
COLNAME_PRIMER1="primerF_seq"
COLNAME_PRIMER2="primerR_seq"

# What proportion of mismatches are you willing to accept when looking for primers?
# Recommended: "0.10"
PRIMER_MISMATCH_PROPORTION="0.10"

################################################################################
# SINGLETONS
################################################################################
# exclude sequences that are found only once? (at the duplicate level)
remove_singletons="YES"  # [YES|NO]


################################################################################
# DEREPLICATION
################################################################################
# Should the sequence ID after dereplication be the output of a hash algorithm?
USE_HASH="NO"

################################################################################
# CLUSTER OTUs
################################################################################
# Would you like to cluster sequences into OTUs based on similarity?
CLUSTER_OTUS="YES"

# What method should be used to cluster OTUs?
cluster_method="swarm" #[ swarm | vsearch | usearch ]

# At what radius of similarity should OTUs be grouped into a cluster?
cluster_radius="1"

# Exclude from the analysis OTUs which are less abundant than what percent?
# Recommendation from Bokulich et al. (2013, Nature Methods): 0.005%
min_OTU_abun="0.005"
# TODO: incorporate into OTU filtering script

################################################################################
# FILTER CHIMERIC SEQUENCES (vsearch)
################################################################################
# Would you like to check for and filter out chimeras?
remove_chimeras="YES"

################################################################################
# TAXONOMIC ANNOTATION
################################################################################
## BLAST ##
PERFORM_BLAST="NO"

# For more information on these parameters, type into a terminal window: blastn -help
# Specify the path to the BLAST database.
# Note this should be a path to any one of three files WITHOUT their extension *.nhr, *.nin, or *.nsq
BLAST_DB='/Users/jimmy.odonnell/NCBI/databases/nt/nt'
# BLAST PARAMETERS
PERCENT_IDENTITY="97"
WORD_SIZE="30"
EVALUE="1e-20"
# number of matches recorded in the alignment:
MAXIMUM_MATCHES="500"
culling_limit="20"


################################################################################
# REANALYSIS
################################################################################
# Would you like to pick up where a previous analysis left off?

# Have the reads already been paired?
ALREADY_PEARED="NO" # YES/NO
PEAR_OUTPUT='/Users/threeprime/Documents/Data/IlluminaData/12S/20140930/Analysis_20141030_2020/1_merged.assembled.fastq.gz'

# Have the merged reads been quality filtered?
ALREADY_FILTERED="NO" # [YES|NO]
FILTERED_OUTPUT='/Users/threeprime/Documents/Data/IlluminaData/12S/20140930/Analysis_20141030_2020/2_filtered_renamed.fasta'


################################################################################
# GENERAL SETTINGS
################################################################################
# Would you like to save every single intermediate file as we go? YES | NO
# recommendation: NO, unless testing or troubleshooting
HOARD="YES"

# Is it ok to rename the sequences within a fasta file?
# This will only remove info about the machine; reads can still be traced back to origin in fastq.
# This will happen after the fastq has been converted to a fasta file at the quality filtering step.
RENAME_READS="NO"

# Make wide format copies of the duplicate and OTU map files?
WIDE_FORMAT="NO"

# Would you like to compress extraneous intermediate files once the analysis is finished? YES/NO
PERFORM_CLEANUP="YES"

# If you want to receive a text message when the pipeline finishes, input your number here:
EMAIL_ADDRESS="4077443377@tmomail.net"
