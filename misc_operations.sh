#!/bin/bash

# To unzip all of the files with .gz extension within a given directory:
my_dir='/path/to/directory'
find "${my_dir}" -type f -name '*.gz' -exec gunzip "{}" \;

# To unzip just the files in the demultiplexed tag folders which had 1 or both tags removed:
DEMULTIPLEXED_DIR='/Users/threeprime/Documents/Data/IlluminaData/12S/20140930/Analysis_20141023_0908/demultiplexed'
find "${DEMULTIPLEXED_DIR}" -type f -name '*prime_tag_rm.fasta.gz' -exec gunzip "{}" \;

# This MIGHT be a way to parallelize that:
# COMPRESSED_FILES=$( find "${DEMULTIPLEXED_DIR}" -type f -name '*prime_tag_rm.fasta.gz' )
# for FILE in "${COMPRESSED_FILES}"; do
# 	( gunzip $FILE ) &
# done
