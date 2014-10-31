#!/bin/bash

# Concatenate demultiplexed reads in order to cluster together

DEMULTIPLEXED_DIR=""
for TAG_SEQ in AACAAC TCACTC TTGAGT; do
  TAG_DIR="${DEMULTIPLEXED_DIR}"/Tag_${TAG_SEQ}
  cat "${TAG_DIR}"/5_primerR_removed.fasta >> ${DEMULTIPLEXED_DIR}/1_newfile.fasta
done



# Concatenate demultiplexed reads in order to cluster together

DEMULTIPLEXED_DIR="/Users/threeprime/Documents/Data/IlluminaData/16S/20141020/Analysis_20141023_1328/demultiplexed"

for (( i=1; i<=11; i++ )); do
  TAG_DIR="${DEMULTIPLEXED_DIR}"/Tag_${i}
  cat "${TAG_DIR}"/both_primer_rem.fasta >> ${DEMULTIPLEXED_DIR}/1_all_sequences.fasta
done
