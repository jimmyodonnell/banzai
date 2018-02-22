#!/usr/bin/env bash

# count the sequences in a fastq or fasta file
# author: Jimmy O'Donnell
# usage: count_seq /path/to/a/file.fastq

count_seq () {
  # detect file type by first character
  CHAR1=$( head -n 1 "${1}" | cut -c-1 )
  if [ "${CHAR1}" == '@' ]; then
    FILETYPE='fastq'
  elif [ "${CHAR1}" == '>' ]; then
    FILETYPE='fasta'
  else
    echo 'unable to properly detect file type'
    return 1
  fi
  
  if [ "${FILETYPE}" == 'fasta' ]; then
    seq_count=$( grep -c '^>' "${1}" )
  else
    fastq_lines=$( cat "${1}" | wc -l )
    seq_count=$(echo "${fastq_lines}" / 4 | bc)
  fi
  echo $seq_count
  return 0
}
