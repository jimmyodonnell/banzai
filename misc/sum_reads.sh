#!/bin/bash

# This script will search all the subdirectories of MY_DIR for any file named TARGET_FILE, and find and return the occurrences of "size=" to the end of the line. It will output the results to a file specified by OUTFILE, which can be modified to show the breakdown of the number of reads per OTU in each file.

MY_DIR='/Users/threeprime/Documents/GoogleDrive/Illumina_Data/12S/run_20140930/Analysis_20141028_1812_JP/demultiplexed'
TARGET_FILE='5_derep.fasta'
OUTFILE='/Users/threeprime/Desktop/outfile.csv'

FILES=$( find "${MY_DIR}" -type f -name "${TARGET_FILE}" )

# write the header line of the csv file
echo "file","sequences" > $OUTFILE

# for each infile, calculate the total number of sequences, and print it along with the shortened file path (everything after the demultiplexed directory)
for INFILE in $FILES; do
	SUM=$( awk -F'=' ' />/ { sub(/;/, "", $NF); print $2 } ' "$INFILE" | awk '{ sum+=$1} END {print sum}' )
	echo "${INFILE#"${MY_DIR}"/}",$SUM >> $OUTFILE
done

# Alt, slower
# awk -F'=' ' />/ { print $2 } ' "$INFILE" | sed 's/;//g' > $OUTFILE

# The old, bad, slow way
# for INFILE in $FILES; do
# 	grep -o -E -H 'size=.*$' "${INFILE}" >> "${OUTFILE}"
# done
# sed 's|'"$MY_DIR"'||g' "${OUTFILE}" > "${OUTFILE%.*}"_tmp.txt
# sed 's|:|,|g' "${OUTFILE%.*}"_tmp.txt > "${OUTFILE}"
# rm "${OUTFILE%.*}"_tmp.txt
