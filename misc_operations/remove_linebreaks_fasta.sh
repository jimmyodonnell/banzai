#!/usr/bin/env bash

# path to input file is given by argument 1
infile="${1}"

# outfile will have same name as input but end with "_l.fasta" (for linear)
outfile="${infile%.*}"_l.fasta

awk '/^>/{print (NR==1)?$0:"\n"$0;next}{printf "%s", $0}END{print ""}' "${infile}" > "${outfile}"
