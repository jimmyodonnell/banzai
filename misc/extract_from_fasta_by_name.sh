#!/usr/bin/env bash

# extract sequences from a fasta file by name
# requires two arguments:
# 1. path to input fasta to draw from
# 2. path to file of names to extract sequences for (one per line)

# sequences must not be wrapped!

echo "Extracting sequences from file:"
infile="${1}"
# infile="/Users/jimmy.odonnell/Desktop/Analysis_20151019_1918/all_lib/OTUs_swarm/OTUs_nochime.fasta"
echo "${infile}"

names="${2}"
# names="/Users/jimmy.odonnell/Desktop/Analysis_20151019_1918/all_lib/OTUs_swarm/OTU_table_filtered_20151105190511/otu_names_samp.txt"
echo "Using the names in file: "
echo "${names}"

echo "Output is:"
outfile="${infile%.*}"_filt.fasta
echo $outfile

grep -f "${names}" "${infile}" -A 1 | sed '/^--$/d' > "${outfile}"

exit
