#!/usr/bin/env bash

primer_F=ACTGGGATTAGATACCCC
primer_R=TAGAACAGGCTCCTCTAG
infile='/Users/threeprime/Desktop/Analysis_20150529_2104/lib1/2_filtered_renamed.fasta'
outfile='cookiecutter_tags.txt'
# grep -E -o '^.{9}ACTGGGATTAGATACCCC' /Users/threeprime/Desktop/Analysis_20150529_2104/lib1/2_filtered.fasta
grep -E -o '^.{6}'"${primer_F}"'' "${infile}" | sed 's/'"${primer_F}"'//' | sort | uniq -c | sort -r > "${outfile}"

