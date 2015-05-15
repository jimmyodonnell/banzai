#!/usr/bin/env bash


################################################################################
# RAW DATA
################################################################################
# What is the file path to the directory containing all of the libraries/reads?
PARENT_DIR="/Users/threeprime/temp_big/test/raw_data/"

# What is the path to the reads?
# READ1='/Users/threeprime/Documents/GoogleDrive/Data_Illumina/16S/run_20150401/libraryA/lib1_R1.fastq.gz'
# READ2='/Users/threeprime/Documents/GoogleDrive/Data_Illumina/16S/run_20150401/libraryA/lib1_R2.fastq.gz'

# This script will generate a directory (folder) containing the output of the script.
# Where do you want this new folder to go?
ANALYSIS_DIRECTORY="/Users/threeprime/Desktop"

# Where is the sequencing pool file? (SEE FORMATTING GUIDELINES IN README!)
SEQUENCING_POOL_DATA="/Users/threeprime/temp_big/run_20150401/20150317_sequencing_pool.csv"

# You can optionally specify a folder into which the script copies a PDF containing some results.
# The pdf is created by default in the analysis folder specified above, but
# if you set this to your DropBox or Google Drive Folder, you can check it out from anywhere.
OUTPUT_PDF_DIR="/Users/threeprime/Desktop"


# TODO grab this from a fragment_size column in the sequencing pool file
# What is the maximum expected length of the fragment of interest, including primers? # AND TAGS?
LENGTH_FRAG="180"

# What is the length of the reads of the Illumina run? (i.e. how long are the sequences in each of the run fastq files (R1 and R2)?)
LENGTH_READ="150"


################################################################################
# MERGE PAIRED READS
################################################################################
# TODO Move PEAR settings here
# For more information on these parameters, type into a terminal window: pear -help


################################################################################
# QUALITY FILTERING
################################################################################
# TODO add quality filtering parameters (mix with PEAR)
# For more information on these parameters, look up the usearch help on Google
# MAXIMUM_EXPECTED_ERRORS
# MINIMUM_SEQUENCE_LENGTH


################################################################################
# HOMOPOLYMERS
################################################################################
# Would you like to remove reads containing runs of consecutive identical bases (homopolymers)?
REMOVE_HOMOPOLYMERS="NO"
# What is the maximum homopolymer length you're willing to accept?
# Reads containing runs of identical bases longer than this will be discarded.
HOMOPOLYMER_MAX="7"


################################################################################
# DEMULTIPLEXING
################################################################################
# Specify the nucleotide sequences that differentiate multiplexed samples ("tags", and in the case of the Kelly Lab, primer tags)
# You can grab these from the file specified above (SEQUENCING_POOL_DATA) by specifying the column name holding tags.
# Or you can specify a text file containing only these tags (choose "NO", and then specify path to the tag file).
# This file should be simply a list of sequences, one per line, of each of the tags, WITH A TRAILING NEWLINE!
# To make a trailing newline, make sure when you open the file, you have hit enter after the final sequence.
READ_TAGS_FROM_SEQUENCING_POOL_DATA="YES" # ["YES"|"NO"] if NO, you must specify the TAG_FILE below
TAG_COLUMN_NAME="tag_sequence"
TAG_FILE='/Users/threeprime/Documents/GoogleDrive/Kelly_Lab_Big/Illumina_Data_Raw/16S/run_20150401/tags_16S.txt'

# How many nucleotides pad the 5' end of the tag sequence?
TAG_Ns="3"
# What is the maximum number of Ns to allow at the end of a sequence before a tag is reached?
# TAG_N_MAX="9" # THIS IS NOT WORKING YET. SET TO DEFAULT 9

# Should demultiplexed samples be concatenated for annotation as a single unit? (Each read can still be mapped back to samples)
# Recommended: YES
CONCATENATE_SAMPLES="YES"

################################################################################
# PRIMER REMOVAL
################################################################################
# Specify the primers used to generate these amplicons.
# As with the multiplex tags, you can grab these from the file SEQUENCING_POOL_DATA.
# In that case, you must indicate the column names of the forward and reverse primers
READ_PRIMERS_FROM_SEQUENCING_POOL_DATA="YES" # ["YES"|"NO"] if NO, you must specify the PRIMER_FILE below
PRIMER_1_COLUMN_NAME="primer_sequence_F"
PRIMER_2_COLUMN_NAME="primer_sequence_R"
# Alternatively, you can specify the path to a fasta file containing the two primers used to generate the amplicons you sequenced:
PRIMER_FILE='/Users/threeprime/Documents/GoogleDrive/Kelly_Lab_Big/Illumina_Data_Raw/16S/run_20150401/primers_16S.fasta'

# What proportion of mismatches are you willing to accept when looking for primers?
# Recommended: "0.10"
PRIMER_MISMATCH_PROPORTION="0.10"



################################################################################
# CLUSTER OTUs
################################################################################
# Would you like to cluster sequences into OTUs based on similarity?
CLUSTER_OTUS="YES"

# What percent similarity must sequences share to be considered the same OTU?
# Note that this must be an integer. Contact me if this is a problem
CLUSTERING_PERCENT="99"


################################################################################
# TAXONOMIC ANNOTATION
################################################################################
## BLAST ##
# For more information on these parameters, type into a terminal window: blastn -help
# Specify the path to the BLAST database.
# Note this should be a path to any one of three files WITHOUT their extension *.nhr, *.nin, or *.nsq
BLAST_DB='/Users/threeprime/temp_big/NCBI/nt_DB/16S_20150511/Metazoa16S20141113.fasta'
# BLAST PARAMETERS
PERCENT_IDENTITY="97"
WORD_SIZE="11"
EVALUE="1e-20"
# number of matches recorded in the alignment:
MAXIMUM_MATCHES="400"

################################################################################
## MEGAN ##
# For more information, see the manual provided with the software
# Specify the path to the MEGAN executable file you want to use.
# Note that in recent versions an executable was not provided; in that case use this:
# /Applications/Megan/MEGAN.app/Contents/MacOS/JavaApplicationStub -g
megan_exec='/Applications/megan/MEGAN'

# What is the lowest taxonomic rank at which MEGAN should group OTUs?
COLLAPSE_RANK1="Family"
MINIMUM_SUPPORT="1"
MINIMUM_COMPLEXITY="0"
TOP_PERCENT="3"
MINIMUM_SUPPORT_PERCENT="0"
MINIMUM_SCORE="140"
LCA_PERCENT="70"
MAX_EXPECTED="1e-25"

# Do you want to perform a secondary MEGAN analysis, collapsing at a different taxonomic level?
PERFORM_SECONDARY_MEGAN="YES"
COLLAPSE_RANK2="Genus"


################################################################################
# REANALYSIS
################################################################################
# Would you like to pick up where a previous analysis left off?
# If reanalyzing existing demultiplexed data, point this variable to the directory storing the individual tag folders.
# EXISTING_DEMULTIPLEXED_DIR='/Users/threeprime/Documents/Data/IlluminaData/16S/20141020/Analysis_20141023_1328/demultiplexed'

# Have the reads already been paired?
ALREADY_PEARED="NO" # YES/NO
PEAR_OUTPUT='/Users/threeprime/Documents/Data/IlluminaData/12S/20140930/Analysis_20141030_2020/1_merged.assembled.fastq.gz'

# Have the merged reads been quality filtered?
ALREADY_FILTERED="NO" # YES/NO
FILTERED_OUTPUT='/Users/threeprime/Documents/Data/IlluminaData/12S/20140930/Analysis_20141030_2020/2_filtered_renamed.fasta'


################################################################################
# GENERAL SETTINGS
################################################################################
# Would you like to compress extraneous intermediate files once the analysis is finished? YES/NO
PERFORM_CLEANUP="YES"

# Is it ok to rename the sequences within a fasta file?
# This will only remove info about the machine; reads can still be traced back to origin in fastq.
# This will happen after the fastq has been converted to a fasta file at the quality filtering step.
RENAME_READS="YES"

# Is the parallel compression utility 'pigz' installed? (Get it here: http://zlib.net/pigz/)
PIGZ_INSTALLED="YES"

# If you want to receive a text message when the pipeline finishes, input your number here:
PHONE_NUMBER="4077443377"
