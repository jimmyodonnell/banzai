#!/usr/bin/env bash

# convert compressed fastq to compressed fasta

# example usage: bash fastqgz_to_fastagz.sh "/Users/threeprime/Desktop/Analysis_20150725_1452/A/1_merged.assembled.fastq.gz"
infile="${1}"

seqtk seq -A "${infile}" | pigz > "${infile%%.*}".fasta.gz
