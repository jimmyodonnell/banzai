#!/bin/bash

# Concatenate demultiplexed reads in order to cluster together

# Option 1:
DEMULTIPLEXED_DIR="/Users/threeprime/Documents/Data/IlluminaData/16S/20141020/Analysis_20141031_0540/demultiplexed"
for FOLDER in $( ls "${DEMULTIPLEXED_DIR}"); do
  DEMULT_FILE="${DEMULTIPLEXED_DIR}"/${FOLDER}/6_primerR_removed.fasta
  cat "${DEMULT_FILE}" >> "${DEMULTIPLEXED_DIR}"/1_demult_concat.fasta
done

# Option 2:
# if the files are zipped:
find "${DEMULTIPLEXED_DIR}" -type f -name '5_primerR_removed.fasta.gz' -exec gunzip "{}" \;
DEMULTIPLEXED_DIR="/Users/threeprime/Documents/Data/IlluminaData/16S/20141020/Analysis_20141031_0540/demultiplexed"
find "${DEMULTIPLEXED_DIR}" -type f -name '5_primerR_removed.fasta' -exec cat "{}" \+ >> "${DEMULTIPLEXED_DIR}"/1_demult_concat.fasta
# recompress
find "${DEMULTIPLEXED_DIR}" -type f -name '5_primerR_removed.fasta' -exec gzip "{}" \;

# rename (clustering will screw up cluster names if there are spaces)
sed 's/M02357:50:000000000-ABLV4:1:1101://g' "${DEMULTIPLEXED_DIR}"/1_demult_concat.fasta > "${DEMULTIPLEXED_DIR}"/1_demult_concat_rename.fasta
sed 's/ 1:N:0:1 /:/g' "${DEMULTIPLEXED_DIR}"/1_demult_concat_rename.fasta > "${DEMULTIPLEXED_DIR}"/1_demult_concat_renamed.fasta
rm "${DEMULTIPLEXED_DIR}"/1_demult_concat_rename.fasta

DEREP_INPUT="${DEMULTIPLEXED_DIR}"/1_demult_concat.fasta
usearch -derep_fulllength "${DEREP_INPUT}" -sizeout -strand both -uc "${DEMULTIPLEXED_DIR}"/derep.uc -output "${DEMULTIPLEXED_DIR}"/7_derep.fasta
DEREP_INPUT="${DEMULTIPLEXED_DIR}"/1_demult_concat_renamed_B.fasta





# Concatenate demultiplexed reads in order to cluster together


for (( i=1; i<=11; i++ )); do
  TAG_DIR="${DEMULTIPLEXED_DIR}"/Tag_${i}
  cat "${TAG_DIR}"/both_primer_rem.fasta >> ${DEMULTIPLEXED_DIR}/1_all_sequences.fasta
done
