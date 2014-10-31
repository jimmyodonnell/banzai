#!/bin/bash

ANALYSIS_DIR='/Users/threeprime/Documents/Data/IlluminaData/12S/20140930/Analysis_20141023_0908/demultiplexed'

for DIR in $(seq 1 1 60)  # ls "$ANALYSIS_DIR"
  do
  FIVEPRIME=$(wc -l < "${ANALYSIS_DIR}"/tag_"${DIR}"/5prime_tag_rm.fasta)
  THREEPRIME=$(wc -l < "${ANALYSIS_DIR}"/tag_"${DIR}"/3prime_tag_rm.fasta)
    echo "${FIVEPRIME}" "${THREEPRIME}" >> tag_chimaeras.txt
  done

# "${ANALYSIS_DIR}"/3_prime_tagged.fasta
# wc -l
