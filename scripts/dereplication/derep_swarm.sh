#!/usr/bin/env bash

# dereplicate with swarm
# usage: bash derep_swarm.sh "/Users/threeprime/Desktop/Analysis_20150727_0004/all_lib/7_no_primers.fasta"
outfile="${1%.*}"_derep_swarm.fasta
size_annotated="${1%.*}"_sizeann.fasta

if grep -q ';size=' "${1}"; then
  echo "Sequences have usearch-style size annotations"
  infile="${1}"
else
  echo "Sequences do not have usearch-style size annotations, and swarm requires them. Adding them..."
  sed '/^>/ s/$/;size=1/' "${1}" > "${size_annotated}"
  infile="${size_annotated}"
fi

n_cores=$(getconf _NPROCESSORS_ONLN)

# swarm --differences 0 --seeds "${outfile}" --threads "${n_cores}" "${infile}"

swarm --differences 0 --seeds "${outfile}" --threads "${n_cores}" --usearch-abundance "${infile}"
