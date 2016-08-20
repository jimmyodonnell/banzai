#!/usr/bin/env bash

# author: Jimmy O'Donnell <jodonnellbio@gmail.com>

# usage: bash correct_NWGC_seqid.sh /path/to/the/directory

# correct the funky sequence IDs in fastq files provided by the UW Northwest Genomics Center (Nickerson Lab)

# Illumina standard format:
# @ <MACHINEID> : <RUNID> : <FLOWCELLID> : <LANE> : <TILE> : <X> : <Y>  (SPACE) <PAIREDREAD> : <FILTERED> : <CONTROL> : <ADAPTERINDEX>

# UW NWGC format:
# @ <flowcell ID> : <lane> : <tile> : <x_pos> : <y_pos> / <read>

the_dir="${1}"

cd "${the_dir}"

out_dir="../renamed"

mkdir "${out_dir}"


the_files=($(find . -type f -name '*.fastq'))

for the_file in "${the_files[@]}"; do

  echo "${the_file}"

  subdir=$(dirname "${the_file##.}")

  subdir_out="${out_dir}""${subdir}"

  mkdir "${subdir_out}"

  seqid_prefix=$(head -n 1 "${the_file}" | grep ^@ | awk -F: '{print $1}' | sort | uniq)
  echo "${seqid_prefix}"

  indexseq=$(basename "${the_file}" | awk -F_ '{print $1}')
  echo "${indexseq}"

  awk -F[:@/] 'BEGIN{OFS=""}; {
    if(/'"${seqid_prefix}"'/)
      print "@machineID:runID:",$2,":",$3,":",$4,":",$5,":",$6," ",$7,":","N",":",0,":","'"${indexseq}"'"
    else
      print $0
  }' "${the_file}" > "${subdir_out}"/"${the_file##*/}"

done
