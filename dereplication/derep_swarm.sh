#!/usr/bin/env bash

# dereplicate with swarm
# usage: bash derep_swarm.sh "/Users/threeprime/Desktop/Analysis_20150727_0004/all_lib/7_no_primers.fasta"
infile="${1}"
outfile="${infile%.*}"_derep_swarm.fasta

n_cores=$(getconf _NPROCESSORS_ONLN)

# swarm --differences 0 --seeds "${outfile}" --threads "${n_cores}" "${infile}"

swarm --differences 0 --seeds "${outfile}" --threads "${n_cores}" --usearch-abundance "${infile}"

