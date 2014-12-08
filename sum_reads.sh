#!/bin/bash

# This script will search all the subdirectories of MY_DIR for any file named TARGET_FILE, and find and return the occurrences of "size=" to the end of the line. It will output the results to a file specified by OUTFILE, which can be modified to show the breakdown of the number of reads per OTU in each file.

MY_DIR='/Users/threeprime/Documents/Data/IlluminaData/12S/run_20140930/Analysis_20141016_1744/demultiplexed'
TARGET_FILE='6_nosingle.fasta'
OUTFILE='/Users/threeprime/Desktop/outfile.txt'

FILES=$( find "${DIRECTORY}" -type f -name "${FILENAME}" )

for INFILE in $FILES; do
	grep -o -E -H 'size=.*$' "${INFILE}" >> "${OUTFILE}"
done

sed 's|'"$DIRECTORY"'||g' "${OUTFILE}" > "${OUTFILE%.*}"_tmp.txt
sed 's|:|,|g' "${OUTFILE%.*}"_tmp.txt > "${OUTFILE}"
rm "${OUTFILE%.*}"_tmp.txt
# awk -F';' ' />/ { print $2 } ' "${INFILE}"
