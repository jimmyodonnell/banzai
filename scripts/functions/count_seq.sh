#!/usr/bin/env bash

# count the sequences in a fastq or fasta file
# author: Jimmy O'Donnell
# usage: count_seq /path/to/a/file.fastq

count_seq () {
  # test for gzip compression
  is_compressed=$( file "${1}" | grep -c 'gzip')

  # detect file type by first character
  if [ "${is_compressed}" = 1 ]; then
    CHAR1=$( gunzip -c "${1}" | head -n 1 | cut -c-1 )
  else
    CHAR1=$( head -n 1 "${1}" | cut -c-1 )
  fi
  
  # set filetype
  if [ "${CHAR1}" == '@' ]; then
    FILETYPE='fastq'
  elif [ "${CHAR1}" == '>' ]; then
    FILETYPE='fasta'
  else
    echo 'unable to properly detect file type'
    return 1
  fi
  
  if [ "${FILETYPE}" == 'fasta' ]; then
    if [ "${is_compressed}" = 1 ]; then
      seq_count=$( gunzip -c "${1}" | grep -c '^>' )
    else
      seq_count=$( grep -c '^>' "${1}" )
    fi
  else
    if [ "${is_compressed}" = 1 ]; then
      fastq_lines=$( gunzip -c "${1}" | wc -l )
    else
      fastq_lines=$( cat "${1}" | wc -l )
    fi
    seq_count=$(echo "${fastq_lines}" / 4 | bc)
  fi
  echo $seq_count
  return 0
}
