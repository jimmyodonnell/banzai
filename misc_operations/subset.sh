#!/usr/bin/env bash

# create a subset of the fastq files stored in a specified directory

# This variable can be changed to grab different numbers of lines (4K lines = 1K sequences)
N_lines=4000

#take argument 1 and set it to variable my_dir
my_dir="${1}"

out_dir="${my_dir}"_sub

# find files with '.fastq.' somewhere in the filename
file_list=($( find "${my_dir}" -type f -name '*.fastq' -o -name '*.fastq.gz' ))

# loop over files found
for current_file in "${file_list[@]}"; do

  echo 'original file is' "${current_file}"

  # get the file path relative to the parent (user-given) directory
  subpath=${current_file##$my_dir}

  # append that relative path, minus the file name, to the output directory
  new_dir="$out_dir${subpath%/*}"

  # make the new directory
  mkdir -p "${new_dir}"

  # get everything after the final / (the filename) of original file
  oldfile="${subpath##*/}"

  # append the "_sub.fastq" to the filename, and add to new directory structure
  newfile="${new_dir}"/"${oldfile%.*}"_sub.fastq
  echo new file is "${newfile}"

  # perform the subsetting:
  # If the extension is .gz, decompress
  if [[ "${current_file}" =~ \.gz$ ]]; then

    gzip -cd "${current_file}" | head -n "${N_lines}" > "${newfile}"

  else

    head -n "${N_lines}" "${current_file}" > "${newfile}"

  fi

  # echo a blank line to space out messages printing to screen
  echo

done
