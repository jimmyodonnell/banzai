#!/usr/bin/env bash


# infile="${1}"
infile="/Users/threeprime/Documents/GoogleDrive/Kelly_Lab/Projects/Lemonade/Data/blast_20151125_1530/blast_results_all.txt"

# output_format="6 qseqid sallseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle"

# query sequence name
query_col=1

# taxonomy column
tax_col=13

# quality column (bitscore or evalue); bitscore is easier to parse in bash
qual_col=12

quer_arr=( $( awk ' { print $'"${query_col}"' } ' "${infile}" | sort | uniq ) )

for iter in "${quer_arr[@]}"; do

    echo "${iter}"
    echo

done

exit


awk '/'"${iter}"'/ print $'""' '
