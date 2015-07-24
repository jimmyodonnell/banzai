#!/usr/bin/env bash

# Run fastqc on any file containing 'fastq' in the filename in a specified directory (and it's subdirectories).

# Example usage: run_fastqc.sh /Users/threeprime/IlluminaRuns/20150717/libraries

my_dir=("$@")
n_cores=$(getconf _NPROCESSORS_ONLN)
output_dir="${my_dir}"/quality_reports

# for output directory at same level as libraries: ${my_dir%/*}/quality_reports



# use this if the fancy bells and whistles fail (thread counts and output directory)
# find "${my_dir}" -type f -name "*.fastq*" | xargs fastqc

find "${my_dir}" -type f -name "*.fastq*" | xargs \
fastqc \
--outdir "${output_dir}" \
--threads "${n_cores}"
