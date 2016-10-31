#!/usr/bin/env bash

# usage:
# source homopolymer_rm.sh <input_fasta_file> <max_homopolymer_length> <outfile>

INPUT_FILE="${1}"
HOMOPOLY_MAX="${2}"
OUTPUT_FILE="${3}"

HomoLineNo="${INPUT_FILE%/*}"/homopolymer_line_numbers.txt
HOMO_FASTA="${INPUT_FILE%/*}"/homopolymeric_reads.fasta

grep -E -i -B 1 -n "(A|T|C|G)\1{$HOMOPOLY_MAX,}" "${INPUT_FILE}" | \
  cut -f1 -d: | \
  cut -f1 -d- | \
  sed '/^$/d' > "${HomoLineNo}"

echo
if [ -s "${HomoLineNo}" ]; then
  awk 'NR==FNR{l[$0];next;} !(FNR in l)' "${HomoLineNo}" "${INPUT_FILE}" > "${OUTPUT_FILE}"
  awk 'NR==FNR{l[$0];next;} (FNR in l)' "${HomoLineNo}" "${INPUT_FILE}" > "${HOMO_FASTA}"
else
  echo "No homopolymers found" > "${NOHOMO_FASTA}"
  mv "${INPUT_FILE}" "${OUTPUT_FILE}"
  echo
fi
