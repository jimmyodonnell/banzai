#!/usr/bin/env bash

# TODO fix subset to not decompress and recompress entire files

# TODO ignore libraries named 'N'

## Annotation/blast
# TODO clean up blast results parsing.
# TODO add automatic output of blast methods.
# TODO consider concatenating all tabular blast hit files into single file


# TODO It appears that the Minimum and Maximum assembly lengths, and the Minimum Overlap lengths being calculated in the banzai.sh script are causing the problem.

# These parameters are calculated from the fragment length and the Read length of the sequence in the fastq file, which is 251.
# Banzai is calculating:
 # Minimum Overlap=176;
 # Minimum assembly length=100;
 # Maximum assembly length=200.
 # I'm using the 150 fragment length that you suggested.
#
# PEAR Minimum Overlap default is 10; Minimum assembly length default is 50; Maximum assembly length default is 0 which disables the restriction.

# @kyamahara: Minimum assembly length = length_fragement - 50;  Maximum assembly length = length_fragment + 50;
#
# OVERLAP_EXPECTED=$(($LENGTH_FRAG - (2 * ($LENGTH_FRAG - $LENGTH_READ) ) ))
# MINOVERLAP=$(( $OVERLAP_EXPECTED / 2 ))
# Where LENGTH_READ = 251

# TODO add sequence table? (this would be a huge file)
possible columns:
original seqid (from sequencer)
modified seqid (redundancies/whitespace removed)
sampleid (lib/tag combo?)
duplicate (duplicate name from dups.fasta)
otu (otu name from otus.fasta)


# TODO add checks for installation of R packages (replace crap in  analyses_prelim.R)

# TODO be aware of potential problems between awk versions on Linux and Mac. On a Mac, the output of `awk --version` is `awk version 20070501` (BSD, I believe)

# TODO add flexibility to stripchart "Reads by sample" in otu analyses

# TODO compress nosingle.txt and 7_no_primers.fasta.derep

# TODO add vsearch clustering

# TODO add library-specific size variable

# TODO An attempt to cause the script to exit if any of the commands returns a non-zero status (i.e. FAILS).

# TODO LIB_TAG_MOD originally contained sort | uniq; this is unnecessary I think

# TODO add `trap "killall background" EXIT` or `trap 'kill $(jobs -p)' EXIT` to kill background processes on exit

# TODO Add decontamination script

# TODO Add normalization

# TODO streamline config file

# TODO Add library names variable (require alpha?)

# TODO write table:
 # - rows: libraries, tags (including rows for whole libraries)
 # - column: numbers for number of sequences: successfully merged, filtered, forward index, both indexes (parallelizable using grep during awk demultiplexing?), primer 1, primer 2, singletons, (dups? otus?)

# TODO Find primer index as a function of 6bp preceeding or following primer seq  `grep -E -o '^.{6}'"${primer_F}"''`

# TODO add --debug / --verbose flag to generate files.
# TODO Potential files to be removed as part of cleanup at the end of script:
# - homopolymer_line_numbers.txt
# - 2_filtered.fasta
# - 1_merged.assembled.fastq
# - 1_merged.assembled_A.fastq
# - 1_merged.assembled_B.fastq



# TODO remove whitespace from sequence labels?
# sed 's/ /_/'


# Copy sequences to fasta files into separate directories based on tag sequence on left side of read
# TODO test for speed against removing the tag while finding it: wrap first tag regex in gsub(/pattern/,""):  awk 'gsub(/^.{0,9}'"$TAG_SEQ"'/,""){if . . .
# TODO MOVE THE VARIABLE ASSIGNMENT TO TOP; MOVE MKDIR TO TOP OF CONCAT IF LOOP
echo $(date +%H:%M) "Concatenating fasta files..."
CONCAT_DIR="$ANALYSIS_DIR"/all_lib
mkdir "${CONCAT_DIR}"
CONCAT_FILE="${CONCAT_DIR}"/1_demult_concat.fasta

# TODO could move this into above loop after demultiplexing?
for CURRENT_LIB in $LIBRARY_DIRECTORIES; do

	LIB_OUTPUT_DIR="${ANALYSIS_DIR}"/${CURRENT_LIB##*/}

	# TODO !!! This will fail if there are underscores in the library names !!!
	# an attempt at making this robust to underscores
	# grep -E -o '_lib_.+?(?=_tag)_tag_.{6}' "${CONCAT_DIR}"/1_demult_concat.fasta | sed 's/_lib_//;s/_tag_/ /' | sort | uniq -c | sort -nr > "${CONCAT_DIR}"/1_demult_concat.fasta.tags

## PRIMER REMOVAL
# TODO: Parallelize cutadapt using gnu parallel: https://github.com/marcelm/cutadapt/issues/157


# TODO rename preliminary to OTU analyses; move analysis script to OTU analysis directory
OUTPUT_PDF="${ANALYSIS_DIR}"/analysis_results_"${START_TIME}".pdf
