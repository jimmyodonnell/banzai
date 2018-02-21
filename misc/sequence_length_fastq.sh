#!/bin/bash

# get length of sequences in a fastq file

INFILE="${1}"

OUTFILE="${1%.*}".seqlen

# write header line
echo 'count length' > "${OUTFILE}"

awk 'NR%4==2 {print length($0)}' "${INFILE}" | sort -nr | uniq -c >> "${OUTFILE}"

# sed + awk to get length of a single read.
# READ_LENGTH=$(sed -n '2p' ${INFILE} | awk '{ print length($0) }')
