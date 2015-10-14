#!/usr/bin/env bash

# This is a shell script to cluster a fasta file of sequences into OTUs using swarm.
# It requires one and only one argument: The absolute path to a fasta file.
# This file should contain only unique sequences - no duplicates.

# example usage: bash cluster_swarm.sh "/Users/threeprime/Desktop/Analysis_20150727_0004/all_lib/no_duplicates.fasta"

# "${1}" is the first argument
infile="${1}"

# maximum number of differences allowed between two amplicons
# two  amplicons will be grouped if they have this many (or fewer) differences
swarm_differences=1

# define output files (these will be in the same directory as the infile)
dir_out="${infile%/*}"/OTUs_swarm
mkdir "${dir_out}"
out_fasta="${dir_out}"/OTUs.fasta
logfile="${dir_out}"/OTUs.log
out_uc="${dir_out}"/OTUs.uc
out_swarm="${dir_out}"/OTUs.swarm
out_stats="${dir_out}"/stats.swarm

# this will automatically find the number of cores on a Unix/Linux computer
n_cores=$(getconf _NPROCESSORS_ONLN)

# execute swarm
swarm \
	--differences "${swarm_differences}" \
	--fastidious \
	--threads "${n_cores}" \
	--output-file "${out_swarm}" \
	--log "${logfile}" \
	--statistics-file "${out_stats}" \
	--uclust-file "${out_uc}" \
	--seeds "${out_fasta}" \
	--usearch-abundance \
	"${infile}" 


exit
