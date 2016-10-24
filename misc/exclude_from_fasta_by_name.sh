#!/usr/bin/env bash

# exclude (remove) sequences from a fasta file by names stored in another file
# requires two arguments:
# 1. path to input fasta to draw from
# 2. path to file of names to extract sequences for (one per line)

# sequences must not be wrapped!

echo "Removing sequences from fasta file:"
infile="${1}"
echo "${infile}"

names="${2}"
echo "Using the names in file: "
echo "${names}"

echo "Output is:"
outfile="${infile%.*}"_filt.fasta
echo $outfile

grep -f "${names}" -v "${infile}" | grep "^>" -A 1 | sed '/^--$/d'  > "${outfile}"

exit
