#!/usr/bin/env bash

# The correct way to do this is to demultiplex, then concatenate, then remove primers and correct orientation for the concatenated sequences.
# This would make more efficient use of cutadapt.
for TAG_SEQ in $TAGS; do
  TAG_DIR="${ANALYSIS_DIR}"/demultiplexed/tag_"${TAG_SEQ}"
  # REMOVE PRIMER SEQUENCES
  # Remove PRIMER1 from the beginning of the reads. NOTE cutadapt1.7+ will accept ambiguities in primers.
  cutadapt -g ^"${PRIMER1_NON}" -e "${PRIMER_MISMATCH_PROPORTION}" --discard-untrimmed "${TAG_DIR}"/4_tagR_removed.fasta > "${TAG_DIR}"/5_primerL1_removed.fasta
  cutadapt -g ^"${PRIMER2_NON}" -e "${PRIMER_MISMATCH_PROPORTION}" --discard-untrimmed "${TAG_DIR}"/4_tagR_removed.fasta > "${TAG_DIR}"/5_primerL2_removed.fasta
  # Remove the reverse complement of PRIMER1 and PRIMER2 from the end of the reads. NOTE cutadapt1.7 will account for anchoring these to the end of the read with $
  cutadapt -a "${PRIMER2RC_NON}" -e "${PRIMER_MISMATCH_PROPORTION}" --discard-untrimmed "${TAG_DIR}"/5_primerL1_removed.fasta > "${TAG_DIR}"/6_primerR1_removed.fasta
  cutadapt -a "${PRIMER1RC_NON}" -e "${PRIMER_MISMATCH_PROPORTION}" --discard-untrimmed "${TAG_DIR}"/5_primerL2_removed.fasta > "${TAG_DIR}"/6_primerR2_removed.fasta
  seqtk seq -r "${TAG_DIR}"/6_primerR2_removed.fasta > "${TAG_DIR}"/6_primerR2_removed_rc.fasta
  cat "${TAG_DIR}"/6_primerR1_removed.fasta "${TAG_DIR}"/6_primerR2_removed_rc.fasta > "${TAG_DIR}"/7_no_primers.fasta
done

DEMULTIPLEXED_DIR="${ANALYSIS_DIR}"/demultiplexed
for FOLDER in $( ls "${DEMULTIPLEXED_DIR}"); do
  DEMULT_FILE="${DEMULTIPLEXED_DIR}"/${FOLDER}/7_no_primers.fasta
  cat "${DEMULT_FILE}" >> "${DEMULTIPLEXED_DIR}"/1_demult_concat.fasta
done
