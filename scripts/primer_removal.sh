#!/usr/bin/env bash

################################################################################
# PRIMER REMOVAL
################################################################################
# requires:
# 1. cutadapt
# 2. seqtk
# 3. revcom function

PRIMER_REMOVAL_INPUT="${1}"
NO_PRIMERS="${2}"
PRIMER1="${3}"
PRIMER2="${4}"
PRIMER_MISMATCH_PROPORTION="${5}"
MINLEN_ROI="${6}"

# intermediate files:
primerL1_removed="${PRIMER_REMOVAL_INPUT%/*}"/5_primerL1_removed.fasta
primerL2_removed="${PRIMER_REMOVAL_INPUT%/*}"/5_primerL2_removed.fasta
primerR1_removed="${PRIMER_REMOVAL_INPUT%/*}"/6_primerR1_removed.fasta
primerR2_removed="${PRIMER_REMOVAL_INPUT%/*}"/6_primerR2_removed.fasta
primerR1_removedRC="${PRIMER_REMOVAL_INPUT%/*}"/6_primerR1_removedRC.fasta

# reverse complement primers
PRIMER1RC=$( revcom "${PRIMER1}" )
PRIMER2RC=$( revcom "${PRIMER2}" )

# count lines in primer removal input
echo $(date +%Y-%m-%d\ %H:%M) "Counting sequences in primer removal input..."
seq_N_demult_concat=$( grep -e '^>' --count "${PRIMER_REMOVAL_INPUT}" )
echo $(date +%Y-%m-%d\ %H:%M) "${seq_N_demult_concat}" "sequences found in primer removal input"
echo

echo $(date +%Y-%m-%d\ %H:%M) "Beginning primer removal..."
# remove primer 1 from left side of sequences
( cutadapt \
	-g ^"${PRIMER1}" \
	-e "${PRIMER_MISMATCH_PROPORTION}" \
	-m "${MINLEN_ROI}" \
	--discard-untrimmed \
	"${PRIMER_REMOVAL_INPUT}" > "${primerL1_removed}" ) &

# remove primer 2 from left side of sequences
( cutadapt \
	-g ^"${PRIMER2}" \
	-e "${PRIMER_MISMATCH_PROPORTION}" \
	-m "${MINLEN_ROI}" \
	--discard-untrimmed \
	"${PRIMER_REMOVAL_INPUT}" > "${primerL2_removed}" ) &

wait

# Remove the reverse complement of primer 1 from the right side of sequences
( cutadapt \
	-a "${PRIMER2RC}"$ \
	-e "${PRIMER_MISMATCH_PROPORTION}" \
	-m "${MINLEN_ROI}" \
	--discard-untrimmed \
	"${primerL1_removed}" > "${primerR1_removed}" ) &

# Remove the reverse complement of primer 2 from the right side of sequences
( cutadapt \
	-a "${PRIMER1RC}"$ \
	-e "${PRIMER_MISMATCH_PROPORTION}" \
	-m "${MINLEN_ROI}" \
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
seqtk seq -r "${primerR1_removed}" > "${primerR1_removedRC}"

# paste together the contents of the files that primers were removed from.
cat "${primerR1_removedRC}" "${primerR2_removed}" > "${NO_PRIMERS}"

# compress or remove input and intermediate files
if [[ "${HOARD}" == "YES" ]]; then
	# compress left primer removal input
	echo $(date +%Y-%m-%d\ %H:%M) "Compressing primer removal files..."
	"${ZIPPER}" "${primerL1_removed}" "${primerL2_removed}" \
		         "${primerR1_removed}" "${primerR2_removed}" "${primerR1_removedRC}"
	echo $(date +%Y-%m-%d\ %H:%M) "Primer removal files compressed."
	echo
else
	rm "${primerL1_removed}" "${primerL2_removed}" \
	   "${primerR1_removed}" "${primerR2_removed}" "${primerR1_removedRC}"
  echo $(date +%Y-%m-%d\ %H:%M) "Primer removal files deleted."
	echo
fi
