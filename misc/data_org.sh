#!/usr/bin/env bash

# move fastq files into subdirectories by their index sequences
# index sequences are the first characters before the '_' (underscore)

# usage: bash data_org.sh /path/to/the/directory

the_dir="${1}"

cd "${the_dir}"

index_seq=($( ls *.fastq | awk -F_ '{print $1}' | sort | uniq ))

for seq in "${index_seq[@]}"; do
  mkdir "${seq}"
  mv "${seq}"*.fastq "${seq}"
done

