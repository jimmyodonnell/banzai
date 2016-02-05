#!/usr/bin/env bash

################################################################################
# EXPECTED ERROR FILTERING (vsearch)
################################################################################
# FILTER READS (This is the last step that uses quality scores, so convert to fasta)

merged_reads_dir="${1}"
# merged_reads_dir="/Users/jimmy.odonnell/Desktop/Ford-PLOSONE/merged_20160201_1408"

Max_Expected_Errors="0.5"

if [ ! -d "${merged_reads_dir}" ]
then
  echo "Argument 1 is not a directory:" "${merged_reads_dir}"
  exit
fi


# get start time
START_TIME=$(date +%Y%m%d_%H%M)
output_dir="${merged_reads_dir%/*}"/filtered_"${START_TIME}"
mkdir "${output_dir}"

# check for multiple cores
n_cores=$(getconf _NPROCESSORS_ONLN)

# Write a log file of output from this script (everything that prints to terminal)
logfile="${output_dir}"/logfile.txt
exec > >(tee "${logfile}") 2>&1

while IFS='' read -r -d '' assembledreads; do
  # echo directory removed is
  # echo "${assembledreads##*/}"
  # echo
  # echo extension removed is
  # echo "${assembledreads%.*}"
  # echo
  nodir="${assembledreads##*/}"
  base="${nodir%.*.*}"
  # echo both removed is
  # echo "${base}"
  # echo vsearch would filter file
  # echo "${assembledreads}"
  # echo
  filtered_output="${output_dir}"/"${base}".fasta
  # echo and make a new file
  # echo "${filtered_output}"
  # echo


  echo $(date +%H:%M) "Filtering merged reads..."
  vsearch \
    --fastq_filter "${assembledreads}" \
    --fastq_maxee "${Max_Expected_Errors}" \
    --fastaout "${filtered_output}" \
    --fasta_width 0

done < <(find "${merged_reads_dir}" -type f -name "*assembled.fastq" -print0)

# clean up the files (compress or delete)
CLEANUP="compress"

if [ "${CLEANUP}" = "compress" ]; then
  echo $(date +%H:%M) "Compressing assembled reads..."
  find "${merged_reads_dir}" -type f -name "*.fastq" -exec pigz "{}" \;
elif [ "${CLEANUP}" = "delete" ]; then
  find "${merged_reads_dir}" -type f -name "*.fastq" -exec rm "{}" \;
fi


echo finished filtering at $(date +%Y%m%d_%H%M)

exit
