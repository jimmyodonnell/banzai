#!/usr/bin/env bash

# make a directory to put all the demultiplexed files in
DEMULTIPLEXED_DIR="${ID1_OUTPUT_DIR}"/demultiplexed
mkdir "${DEMULTIPLEXED_DIR}"

# Copy sequences to fasta files into separate directories based on index sequence on left side of read
# TODO test for speed against removing the index while finding it: wrap first index regex in gsub(/pattern/,""):  awk 'gsub(/^.{0,9}'"$IND_SEQ"'/,""){if . . .
# 20150522 changed {0,9} to {3} to eliminate flexibility (that could result in a read being assigned to >1 sample)
# awk '/^.{0,9}'"$IND_SEQ"'/{if (a && a !~ /^.{0,9}'"$IND_SEQ"'/) print a; print} {a=$0}' "${DEMULTIPLEX_INPUT}" > "${ID2_DIR}"/1_tagL_present.fasta ) &

echo $(date +%Y-%m-%d\ %H:%M) "Demultiplexing: removing secondary index sequences and adding to sequence ID in library" "${CURRENT_LIB##*/}""..."
for IND_SEQ in "${ID2S[@]}"; do
(	ID2_DIR="${DEMULTIPLEXED_DIR}"/tag_"${IND_SEQ}"
  mkdir "${ID2_DIR}"
  demult_file_L="${ID2_DIR}"/1_tagL_removed.fasta
  demult_file_R="${ID2_DIR}"/2_notags.fasta

  # Left side secondary index
  awk 'gsub(/^.{3}'"$IND_SEQ"'/,"") {
    if (a && a !~ /^.{3}'"$IND_SEQ"'/)
      print a;
    print
  } {a=$0}' "${DEMULTIPLEX_INPUT}" > "${demult_file_L}"

  # Right side secondary index
  ID2_RC=$( echo ${IND_SEQ} | tr "[ATGCatgcNn]" "[TACGtacgNn]" | rev )
  awk 'gsub(/'"$ID2_RC"'.{3}$/,"") {
    if (a && a !~ /'"$ID2_RC"'.{3}$/)
      print a ";ID2=""'"$IND_SEQ"'";
    print
  } {a = $0}' "${demult_file_L}" > "${demult_file_R}"

  echo "${CURRENT_ID1_NAME}" "${IND_SEQ}" $(wc -l "${demult_file_L}" | \
    awk '{ print ($1/2) }') $(wc -l "${demult_file_R}" | \
    awk '{ print ($1/2)}') >> "${INDEX_COUNT}" ) &

done

wait
