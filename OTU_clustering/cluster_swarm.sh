#!/usr/bin/env bash

# This is a shell script to cluster a fasta file of sequences into OTUs using swarm.
# It requires one and only one argument: The absolute path to a fasta file.
# This file should contain only unique sequences - no duplicates.

# example usage: bash cluster_swarm.sh "/Users/threeprime/Desktop/Analysis_20150727_0004/all_lib/no_duplicates.fasta"

# "${1}" is the first argument
infile="${1}"

# maximum number of differences allowed between two amplicons
# two  amplicons will be grouped if they have this many (or fewer) differences
cluster_radius="$(( 100 - ${CLUSTERING_PERCENT} ))"
# swarm_differences=1

# define output files (these will be in the same directory as the infile)
OTU_dir="${infile%/*}"/OTUs_swarm
mkdir "${OTU_dir}"
OTU_fasta="${OTU_dir}"/OTUs.fasta
logfile="${OTU_dir}"/OTUs.log
out_swarm="${OTU_dir}"/OTUs.swarm
out_stats="${OTU_dir}"/stats.swarm
dup_otu_map="${OTU_dir}"/dups_to_otus.csv
BLAST_INPUT="${OTU_fasta}"
OTU_table="${OTU_dir}"/OTU_table.csv

# this will automatically find the number of cores on a Unix/Linux computer
n_cores=$(getconf _NPROCESSORS_ONLN)

# execute swarm
swarm \
	--differences "${cluster_radius}" \
	--fastidious \
	--threads "${n_cores}" \
	--output-file "${out_swarm}" \
	--log "${logfile}" \
	--statistics-file "${out_stats}" \
	--seeds "${OTU_fasta}" \
	--usearch-abundance \
	"${infile}"

	# Note from swarm author: The -u option also slows does swarm a lot, do not use it if you don't really need results in that format.
	# out_uc="${dir_out}"/OTUs.uc
	# --uclust-file "${out_uc}" \



# make a file mapping each duplicate to its OTU so OTU contingency table can be generated.
# infile should be the output file from swarm -- currently called "OTUs.swarm"

awk 'BEGIN{
            print "Query,Match"
          }
          {
            c = split($0, s);
            for(n=1; n<=c; ++n)
            print s[n] "," $1
          }' "${out_swarm}" |\
sed 's/;size=[0-9]*;//g' > "${dup_otu_map}"

# exit
