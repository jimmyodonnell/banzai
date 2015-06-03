#!/usr/bin/env bash

# usage: bash count_sequences.sh '/Path/to/directory' > outputfile.txt

args=("$@")

find "${args[0]}" -type f -name '*.fasta.gz' -exec gunzip "{}" \;

# files=$(find "${args[0]}" -type f -name '*.fasta.gz' )
