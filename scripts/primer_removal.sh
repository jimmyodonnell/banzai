#!/usr/bin/env bash

################################################################################
# PRIMER REMOVAL
################################################################################
# count lines in primer removal input
echo $(date +%Y-%m-%d\ %H:%M) "Counting sequences in primer removal input..."
seq_N_demult_concat=$( grep -e '^>' --count "${CONCAT_FILE}" )
echo $(date +%Y-%m-%d\ %H:%M) "${seq_N_demult_concat}" "sequences found in primer removal input" #"${CONCAT_FILE}"
echo

echo $(date +%Y-%m-%d\ %H:%M) "Beginning primer removal..."
# remove primer 1 from left side of sequences
primerL1_removed="${CONCAT_DIR}"/5_primerL1_removed.fasta
( cutadapt \
	-g ^"${PRIMER1}" \
	-e "${PRIMER_MISMATCH_PROPORTION}" \
	-m "${LENGTH_ROI_HALF}" \
	--discard-untrimmed \
	"${CONCAT_FILE}" > "${primerL1_removed}" ) &

# remove primer 2 from left side of sequences
primerL2_removed="${CONCAT_DIR}"/5_primerL2_removed.fasta
( cutadapt \
	-g ^"${PRIMER2}" \
	-e "${PRIMER_MISMATCH_PROPORTION}" \
	-m "${LENGTH_ROI_HALF}" \
	--discard-untrimmed \
	"${CONCAT_FILE}" > "${primerL2_removed}" ) &

wait

# compress left primer removal input
echo $(date +%Y-%m-%d\ %H:%M) "Compressing left primer removal input..."
"${ZIPPER}" "${CONCAT_DIR}"/1_demult_concat.fasta
echo $(date +%Y-%m-%d\ %H:%M) "Left primer removal input compressed."
echo

# check for cutadapt/primer removal success.
if [[ ! -s "${primerL1_removed}" ]]; then
  echo 'ERROR: cutadapt did not process reads correctly. This file is empty or absent:'
	echo "${primerL1_removed}"
  echo 'Aborting script'
  exit
fi
# check for cutadapt/primer removal success.
if [[ ! -s "${primerL2_removed}" ]]; then
  echo 'ERROR: cutadapt did not process reads correctly. This file is empty or absent:'
	echo "${primerL2_removed}"
  echo 'Aborting script'
  exit
fi

# Remove the reverse complement of primer 1 from the right side of sequences
primerR1_removed="${CONCAT_DIR}"/6_primerR1_removed.fasta
( cutadapt \
	-a "${PRIMER2RC}"$ \
	-e "${PRIMER_MISMATCH_PROPORTION}" \
	-m "${LENGTH_ROI_HALF}" \
	--discard-untrimmed \
	"${primerL1_removed}" > "${primerR1_removed}" ) &

# Remove the reverse complement of primer 2 from the right side of sequences
primerR2_removed="${CONCAT_DIR}"/6_primerR2_removed.fasta
( cutadapt \
	-a "${PRIMER1RC}"$ \
	-e "${PRIMER_MISMATCH_PROPORTION}" \
	-m "${LENGTH_ROI_HALF}" \
	--discard-untrimmed \
	"${primerL2_removed}" > "${primerR2_removed}" ) &

wait

# check for cutadapt/primer removal success.
if [[ ! -s "${primerR1_removed}" ]]; then
	echo 'ERROR: cutadapt did not process reads correctly. This file is empty or absent:'
	echo "${primerR1_removed}"
	echo 'Aborting script'
	exit
fi
# check for cutadapt/primer removal success.
if [[ ! -s "${primerR2_removed}" ]]; then
	echo 'ERROR: cutadapt did not process reads correctly. This file is empty or absent:'
	echo "${primerR2_removed}"
	echo 'Aborting script'
	exit
fi

echo

# Reverse-complement the sequences in which the RC of primer 1 was found on the right side
echo $(date +%Y-%m-%d\ %H:%M) "Correcting sequence orientation..."
seqtk seq -r "${CONCAT_DIR}"/6_primerR1_removed.fasta > "${CONCAT_DIR}"/6_primerR1_removedRC.fasta

# paste together the contents of the files that primers were removed from.
DEREP_INPUT="${CONCAT_DIR}"/7_no_primers.fasta
cat "${CONCAT_DIR}"/6_primerR1_removedRC.fasta "${CONCAT_DIR}"/6_primerR2_removed.fasta > "${DEREP_INPUT}"

# check that it worked (derep input / no primers)
if [[ ! -s "${DEREP_INPUT}" ]] ; then
    echo 'ERROR: Input file for dereplication is empty or absent.'
    echo 'This will cause problems for all remaining steps, so script will exit.'
    exit
fi

echo
