#!/bin/bash

############################################################################
# COUNT OCCURRENCES PER SAMPLE (LIBRARY + TAG) PER DUPLICATE
## Filename was "write_otu_table.sh"

INFILE='/Users/threeprime/Documents/Data/IlluminaData/16S/run_20141113_time_series/all_libraries/nosingle'
PRIMER_TAGS_LIB='/Users/threeprime/Documents/Data/IlluminaData/16S/run_20141113_time_series/tags_16S_lib.txt'
TAGS=$(tr '\n' ' ' < "${PRIMER_TAGS_LIB}" )

for TAG_SEQ in $TAGS; do
	( awk 'BEGIN {print "'$TAG_SEQ'" ; FS ="_tag_'${TAG_SEQ}'" } { print NF -1 }' "${INFILE}" > "${INFILE%/*}"/"${TAG_SEQ}".dup ) &
done
wait

# Write a csv file of the number of occurrences of each duplicate sequence per tag.
find "${DEREP_INPUT%/*}" -type f -name '*.dup' -exec paste -d, {} \+ | awk '{ print "DUP_" NR-1 "," $0 }' > "${DEREP_INPUT%/*}"/dups.csv
rm ${DEREP_INPUT%/*}/*.dup

# to just print the IDs of the sequences in a FASTA file
# awk ' NR % 2 ' "${FASTA_FILE}" > "${FASTA_FILE%/*}"/seqIDs.txt
