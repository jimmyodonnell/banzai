#!/usr/bin/env bash

# create a subset of the fastq files stored in a specified directory

# This variable can be changed to grab different numbers of lines (4K lines = 1K sequences)
N_lines=4000

#take argument 1 and set it to variable my_dir
my_dir="${1}"

out_dir="${my_dir}"_sub

# find files with '.fastq.' somewhere in the filename
file_list=($( find "${my_dir}" -type f -name '*.fastq' -o -name '*.fastq.gz' ))

# test whether pigz is installed
if command -v pigz >/dev/null 2>&1; then
  echo "pigz is installed"
  zipper="pigz"
else
  echo "pigz not installed"
  zipper="gzip"
fi

echo

# loop over files found
for current_file in "${file_list[@]}"; do

  # If the extension is .gz, decompress
  if [[ "${current_file}" =~ \.gz$ ]]; then

    echo 'decompressing' "${current_file}"
    "${zipper}" -d "${current_file}"

    my_fastq="${current_file%.gz}"

  else

    my_fastq="${current_file}"

  fi

  echo 'fastq file is' "${my_fastq}"

  # get the file path relative to the parent (user-given) directory
  subpath=${my_fastq##$my_dir}

  # append that relative path, minus the file name, to the output directory
  new_dir="$out_dir${subpath%/*}"

  # make the new directory
  mkdir -p "${new_dir}"

  # get everything after the final / (the filename) of original file
  oldfile="${subpath##*/}"

  # append the "_sub.fastq" to the filename, and add to new directory sturcture
  newfile="${new_dir}"/"${oldfile%.*}"_sub.fastq
  echo new file is "${newfile}"

  head -n "${N_lines}" "${my_fastq}" > "${newfile}"

  # if the input file was compressed, compress it again.
  if [[ "${current_file}" =~ \.gz$ ]]; then

    "${zipper}" "${my_fastq}"

  fi

  # echo a blank line to space out messages printing to screen
  echo

done
