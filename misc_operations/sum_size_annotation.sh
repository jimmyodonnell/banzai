#!/usr/bin/env bash

# usage: bash sum_size_annotation.sh /path/to/file.fasta

# sum sequence size:
awk -F';|=' '/^>/ { sum+=$3 } END { print sum} ' "${1}"
