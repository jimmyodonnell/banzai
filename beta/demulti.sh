#!/usr/bin/env bash

# Author: Jimmy O'Donnell jodonnellbio@gmail.com

# Usage    bash demulti.sh -i infile.fasta -s 4 -l 6

# This script will:
# 1. read a fasta file (sequences must NOT occupy multiple lines!)
# 2. trim sequences at specified location from both sides (left and right)
# 3. place the sequences in the fasta sequence header
# 4. print results to STDOUT (aka the screen; you can redirect/pipe from there)

# It requires you to provide:
# -i | --infile : a path to a fasta file
# -s | --start  : at which position does the index start? (1 or more)
# -l | --length : how many bases long is the index? (1 or more)

# and optionally, you can indicate
# --leftonly : indexes are only on left side

# TODO: check for line wraps in sequences

while [[ $# > 1 ]]
do
key="$1"

case $key in
    -i|--infile)
    INFILE="$2"
    shift # past argument
    ;;
    -s|--start)
    INDEXSTART="$2"
    shift # past argument
    ;;
    -l|--length)
    INDEXLENGTH="$2"
    shift # past argument
    ;;
    --leftonly)
    LEFTONLY=YES
    ;;
    *)
            # unknown option
    ;;
esac
shift # past argument or value
done


SEQ_START=$(( INDEXSTART + INDEXLENGTH ))
INDEX_BOTH=$(( (2*${SEQ_START}) - 2 ))

# echo INPUT FILE           = "${INFILE}"
# echo INDEX START POSITION = "${INDEXSTART}"
# echo INDEX LENGTH         = "${INDEXLENGTH}"
# echo $SEQ_START
# echo $INDEX_BOTH

if [[ -n $1 ]]; then
	awk 'BEGIN {OFS=""};
	  /^[ACTG]/ {
	  print a \
		";ID2A=" substr($0,'${INDEXSTART}','${INDEXLENGTH}')\
		";ID2B=NONE";
		print substr($0, '${SEQ_START}', length($0))
	  } {a=$0}' "${INFILE}"
else
	awk 'BEGIN {OFS=""};
	  /^[ACTG]/ {
	  print a \
		";ID2A=" substr($0,'${INDEXSTART}','${INDEXLENGTH}')\
		";ID2B=" substr($0,length($0)+2-'${SEQ_START}','${INDEXLENGTH}');
		print substr($0, '${SEQ_START}', length($0)-'${INDEX_BOTH}')
	  } {a=$0}' "${INFILE}"
fi
