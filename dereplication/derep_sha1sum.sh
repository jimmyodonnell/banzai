#!/usr/bin/env bash

# dereplication

# usage: bash derep.sh '/Path/to/file.fasta'

infile="${1}"
outfile="${infile%.*}"_derep_bash.fasta

grep -v "^>" "${1}" | \
grep -v [^ACGTacgt] | sort -d | uniq -c | \
while read abundance sequence ; do
    hash=$(printf "${sequence}" | sha1sum)
    hash=${hash:0:40}
    printf ">%s_%d_%s\n" "${hash}" "${abundance}" "${sequence}"
done | sort -t "_" -k2,2nr -k1.2,1d | \
sed -e 's/\_/\n/2' > "${outfile}"

