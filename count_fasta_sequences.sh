#!/bin/bash

# This script will find any file named TARGET_FILE in any subdirectory of directory MY_DIR, and count the number of lines containing '>' (i.e. the number of sequences in a FASTA file).

MY_DIR='/Users/threeprime/Documents/Data/IlluminaData/12S/run_20140930/Analysis_20141016_1744/demultiplexed'
TARGET_FILE='3prime_tag_rm.fasta'
OUTPUT_FILE='/Users/threeprime/Desktop/output.txt'

FILES=$( find "${MY_DIR}" -type f -name "${TARGET_FILE}" )

for FILE in $FILES; do
	grep -c -H '>' "${FILE}" >> "${OUTPUT_FILE}"
done
