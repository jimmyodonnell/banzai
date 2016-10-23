bas#!/usr/bin/env bash

# output from illumina machines is fastq files.
# For paired-end sequencing runs, there will be two files output that must later be merged.

# suggested usage:
# bash merge_files_in_one_dir.sh "path/to/folder" "R1" "R2"

# this script requires three arguments:

# 1. a directory holding some fastq files

# 2. a string specifying first read identifier
readidentifier1="${2}"
# readidentifier1="R1"

# 3. a string specifying second read identifier
readidentifier2="${3}"
# readidentifier2="R2"

# read and check argument 1
mydir="${1}"
# mydir="/Users/jimmy.odonnell/Desktop/Ford-PLOSONE/data_uncomp"

if [ ! -d "${mydir}" ]
then
  echo "Argument 1 is not a directory:" "${mydir}"
  exit
fi


# 4. (optional) output directory
# by default, files will be placed in a new directory at the same level as the input
if [ -z "${4}" ]; then
  output_parent="${mydir%/*}"
else
  output_parent="${4}"
fi

# get start time
START_TIME=$(date +%Y%m%d_%H%M)
output_dir="${output_parent}"/merged_"${START_TIME}"
mkdir "${output_dir}"


# Write a log file of output from this script (everything that prints to terminal)
logfile="${output_dir}"/logfile.txt
exec > >(tee "${logfile}") 2>&1

# check for multiple cores
n_cores=$(getconf _NPROCESSORS_ONLN)

while IFS='' read -r -d '' firstread; do
  # echo first read is
  # echo "${firstread}"
  # echo second read is
  basename="${firstread##*/}"
  secondread="${firstread%/*}"/${basename/$readidentifier1/$readidentifier2}
  # echo "${secondread}"

  pear \
    --forward-fastq "${firstread}" \
    --reverse-fastq "${secondread}" \
    --output "${output_dir}"/"${basename%.*}" \
    -j $n_cores

done < <(find "${mydir}" -type f -name "*${readidentifier1}*.fastq" -print0)

echo finished merging at $(date +%Y%m%d_%H%M)

exit
