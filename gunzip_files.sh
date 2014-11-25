#!/bin/bash

# This is a way to gunzip all files matching certain criteria within a given directory:

# Enter your search criteria after the argument ' -name '.
# For example:
# -name '*.fasta.gz' # to find any gzipped fasta file
# -name '*.fastq.gz' # to find any gzipped fastq file
# -name 'both_primer_rem.fasta.gz' # to find any file named both_primer_rem.fasta.gz

MY_DIR='/Users/path/to/the/directory'
find "${MY_DIR}" -type f -name '*.fasta.gz' -exec gunzip "{}" \;
