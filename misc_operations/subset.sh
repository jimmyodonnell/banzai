#!/usr/bin/env bash

# create a subset of the fastq files stored in a specified directory

# This variable can be changed to grab different numbers of lines (4K lines = 1K sequences)
N_lines=4000

#take argument 1 and set it to variable my_dir
my_dir="${1}"

# find files with '.fastq.' somewhere in the filename
file_list=($( find "${my_dir}" -type f -name '*.fastq*' ))

# test whether pigz is installed
which pigz

if [ $? -eq 0 ]; then
  echo "pigz is installed"
  zipper="pigz"
else
  echo "pigz not installed"
  zipper="gzip"
fi

# loop over files found
for file in "${file_list[@]}"; do

  # If the extension is .gz
  if [[ "${file}" =~ \.gz$ ]]; then

    "${zipper}" -d "${file}"

    my_fastq="${file%.*}"

  else

    my_fastq="${file}"

  fi

  head -n "${N_lines}" "${my_fastq}" > "${my_fastq%.*}"_1K.fastq

  # if the input file was compressed, compress it again.
  if [[ "${file}" =~ \.gz$ ]]; then

    "${zipper}" "${my_fastq}"

  fi

done
