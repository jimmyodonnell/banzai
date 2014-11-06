#!/bin/bash

DEMULTIPLEXED_DIR='/Users/threeprime/Documents/Data/IlluminaData/12S/20140930/Analysis_20141023_0908/demultiplexed'
echo "five_prime_tagged" "three_prime_tagged">> "${DEMULTIPLEXED_DIR}"/tag_chimaeras.txt
for TAG_DIR in $( ls "${DEMULTIPLEXED_DIR}" ); do
  FIVEPRIME=$(wc -l < "${DEMULTIPLEXED_DIR}"/"${TAG_DIR}"/5prime_tag_rm.fasta)
  THREEPRIME=$(wc -l < "${DEMULTIPLEXED_DIR}"/"${TAG_DIR}"/3prime_tag_rm.fasta)
  echo "${FIVEPRIME}" "${THREEPRIME}" >> "${DEMULTIPLEXED_DIR}"/tag_chimaeras.txt
done
