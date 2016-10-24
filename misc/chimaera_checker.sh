#!/usr/bin/env bash

DEMULTIPLEXED_DIR='/Users/threeprime/Documents/Data/IlluminaData/12S/20140930/Analysis_20141023_0908/demultiplexed'
echo "tag_directory" "five_prime_tagged" "three_prime_tagged">> "${DEMULTIPLEXED_DIR}"/tag_chimaeras.txt
for TAG_DIR in $( ls "${DEMULTIPLEXED_DIR}" ); do
  FIVEPRIME=$(wc -l < "${DEMULTIPLEXED_DIR}"/"${TAG_DIR}"/5prime_tag_rm.fasta)
  THREEPRIME=$(wc -l < "${DEMULTIPLEXED_DIR}"/"${TAG_DIR}"/3prime_tag_rm.fasta)
  echo "${TAG_DIR}" "${FIVEPRIME}" "${THREEPRIME}" >> "${DEMULTIPLEXED_DIR}"/tag_chimaeras.txt
done

# sorting (imperfect)
# (head -n 1 <"${DEMULTIPLEXED_DIR}"/tag_chimaeras.txt> && tail -n +3 <"${DEMULTIPLEXED_DIR}"/tag_chimaeras.txt> | sort) > newfile
