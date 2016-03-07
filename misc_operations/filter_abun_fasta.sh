#!/usr/bin/env bash

# this script reads a fasta file, and outputs a fasta file excluding sequences whose abundance is less than a specified threshold
# takes two arguments:
# 1. the path to a fasta file
# 2. the minimum abundance required to keep a sequence in output
# example usage:
# bash "/path/to/this/script/filter_abun_fasta.sh" "/path/to/a/file.fasta" "100"

# assign argument 1 to a variable
infile="${1}"

if [ ! -e "${infile}" ]; then
  echo "Error: Specified infile does not exist:"
  echo "${infile}"
  exit 1
fi

# assign argument 2 to a variable
min_abun="${2}"

outfile=${infile%.*}_filt_abun.fasta
# outfile_toofew=${infile%.*}_filt_abun_toofew.fasta


# check if min_abun is indeed an integer
case $min_abun in
    ''|*[!0-9]*)
      echo "Error: minimum abundance argument ('""${min_abun}""') is not an integer" >&2
      exit 1
    ;;

    *) echo "Removing sequences with fewer than ""${min_abun}"" occurrences from specified file:"
    echo "${infile}"

    ;;

esac

# TODO print sequences with too few occurances to separate file
# see: http://stackoverflow.com/questions/13077631/is-it-possible-to-print-different-lines-to-different-output-files-using-awk
awk -v min="${min_abun}" '\
  BEGIN {
    FS = ";size="
  } \
  { if( $2 >= min )\
    { print;getline;print } \
  }' "${infile}" > "${outfile}"

echo
echo "Sequences occurring ""${min_abun}"" or more times written to file:"
echo "${outfile}"
# echo
# echo "Sequences occurring fewer than ""${min_abun}"" times written to file:"
# echo "${outfile_toofew}"
