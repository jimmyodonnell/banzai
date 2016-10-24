#!/usr/bin/env bash

# Demultiplexed sequencing data often comes back from a sequencing facility buried in a deep file hierarchy. For example:
#  /PE2x151_RK-P_RyanKelly032615-21307490/1-23491099/Data/Intensities/BaseCalls/library1_S1_L001_R1_001.fastq.gz

# It is useful to move those files up in the file hierarchy, so that all fastq files are at level 1 instead of level 5 (as above)


# As written (mv "${dir}"/*/*/*/*.fastq.gz "${dir}"), it will move the files nested 5 levels below the parent directory ("level0") into folders at "level1".
# /level0/level1/level2/level3/level4/file.fastq.gz
# /level0/level1/file.fastq.gz

# YOU SHOULD TEST THIS SCRIPT FIRST ON A TEMPORARY DIRECTORY

# usage:
# bash "path/to/organize_files.sh" "/path/to/directory/with/subdirectories"

parent_dir="${1}"

# find "${parent_dir}" -type d -depth 1

for dir in "${parent_dir}"/*; do
	mv "${dir}"/*/*/*/*.fastq.gz "${dir}"
# 	rm -rf "${dir}"/Data
done

exit


