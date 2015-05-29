#!/usr/bin/env bash

# usage: bash count_sequences.sh '/Path/to/directory' > outputfile.txt
args=("$@")

start_time=$(date +%Y%m%d_%H%M%S)

outfile=$( echo "${args[0]}"/"counted_sequences_""${start_time}"".txt" )

files=$(find "${args[0]}" -type f -name *.fasta )
declare -a files_array=($files)
for i in "${files_array[@]}"; do
	seqs=$( grep -c '^>' "${i}" )
	echo ${i/"${args[0]}"/} $seqs | sed 's|[/ ]|,|g' | sed 's|^,||' >> "${outfile}"
done

# echo ${i/"$args"/}


##################################
# Strike all of this below.
# this script will count the number of sequences in each of the subdirectories of a folder ("DEMULTIPLEXED_DIR")
# It does this by counting the number of lines containing a '>' character
# The following files will be read from each folder:
# 5prime_tag_rm.fasta
# 3prime_tag_rm.fasta
# 5prime_primer_rm.fasta
# both_primer_rem.fasta
# 5_derep.fasta
# 6_nosingle.fasta
# 7_OTUs.fasta

# Specify the directory to search.
# DEMULTIPLEXED_DIR='/Users/threeprime/Documents/GoogleDrive/Illumina_Data/12S/run_20140930/Analysis_20141028_1812_JP/demultiplexed'
# 
# # Specify the output file location
# OUTPUT_FILE="${DEMULTIPLEXED_DIR}"/sequence_counts.csv
# 
# echo "tag_directory","five_prime_tagged","three_prime_tagged","five_prime_primed","three_prime_primed","derep","no_singletons","OTUs">> "${OUTPUT_FILE}"
# for TAG_DIR in "${DEMULTIPLEXED_DIR}"/*/; do
#   TAG_FIVEPRIME=$(grep -c '>' "${TAG_DIR}"/5prime_tag_rm.fasta)
#   TAG_THREEPRIME=$(grep -c '>' "${TAG_DIR}"/3prime_tag_rm.fasta)
#   PRIMER_FIVEPRIME=$(grep -c '>' "${TAG_DIR}"/5prime_primer_rm.fasta)
#   PRIMER_THREEPRIME=$(grep -c '>' "${TAG_DIR}"/both_primer_rem.fasta)
#   DEREP=$(grep -c '>' "${TAG_DIR}"/5_derep.fasta)
#   NOSING=$(grep -c '>' "${TAG_DIR}"/6_nosingle.fasta)
#   OTUS=$(grep -c '>' "${TAG_DIR}"/7_OTUs.fasta)
#   DIR_PATH="${TAG_DIR%/*}"
#   echo "${DIR_PATH##*/}","${TAG_FIVEPRIME}","${TAG_THREEPRIME}","${PRIMER_FIVEPRIME}","${PRIMER_THREEPRIME}","${DEREP}","${NOSING}","${OTUS}">> "${OUTPUT_FILE}"
# done

# sorting (imperfect)
# (head -n 1 <"${DEMULTIPLEXED_DIR}"/tag_chimaeras.txt> && tail -n +3 <"${DEMULTIPLEXED_DIR}"/tag_chimaeras.txt> | sort) > newfile
