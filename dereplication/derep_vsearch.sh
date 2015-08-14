#!/usr/bin/env bash

# dereplicate with vsearch

# specify input file
infile="${1}"

# perform parameter expansion to manipulate output file names
out_dir="${infile%/*}"/derep_vsearch
mkdir "${out_dir}"

outfile="${out_dir}"/derep.fasta
outfile_uc="${out_dir}"/derep.uc
logfile="${out_dir}"/derep.log

# execute vsearch
vsearch --derep_fulllength "${infile}" --fasta_width 0 --sizeout --uc "${outfile_uc}" --log "${logfile}" --output "${outfile}"
