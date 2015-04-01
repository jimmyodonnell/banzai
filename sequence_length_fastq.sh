#!/bin/bash

# get length of sequences in a fastq file

OUTFILE='/Users/threeprime/Desktop/outfile3.txt'
INFILE='/Users/threeprime/Documents/GoogleDrive/Data_Illumina/16S/run_20141113_time_series/library3/Analysis_20141113_2055/1_merged.assembled.fastq'
awk 'NR%4==2 {print length($0)}' "${INFILE}" > "${OUTFILE}"

# sed + awk to get length of a single read.
READ_LENGTH=$(sed -n '2p' ${INFILE} | awk '{ print length($0) }')
