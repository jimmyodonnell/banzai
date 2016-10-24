#!/usr/bin/env bash

# this script reads a fasta file, and outputs a fasta file excluding sequences shorter than the specified integer
# takes two arguments:
# 1. the path to a fasta file
# 2. the minimum length of sequences to be output
# example usage:
# bash "/path/to/this/script/filter_length_fasta.sh" "/path/to/a/file.fasta" "100"

# assign argument 1 to a variable
infile="${1}"

if [ ! -e "${infile}" ]; then
  echo "Error: Specified infile does not exist:"
  echo "${infile}"
  exit 1
fi

outfile=${infile%.*}_filt_len.fasta

# assign argument 2 to a variable
min_length="${2}"

# check if min_length is indeed an integer
case $min_length in
    ''|*[!0-9]*)
      echo "Error: minimum length argument ('""${min_length}""') is not an integer" >&2
      exit 1
    ;;

    *) echo "Removing sequences shorter than ""${min_length}"" bp from specified file:"
    echo "${infile}"

    ;;

esac

awk -v min="${min_length}" '\
  BEGIN {
    RS = ">" ; ORS = ""
  } \
  length($2) >= min { \
    print ">"$0 \
  }' "${infile}" > "${outfile}"

echo "Sequences longer than ""${min_length}"" bp written to file:"
echo "${outfile}"
