#!/usr/bin/env bash

if [ "${ALREADY_FILTERED}" = "YES" ]; then
  echo "Using existing filtered reads in file $FILTERED_OUTPUT"
else
  FILTERED_OUTPUT="${LIB_OUTPUT_DIR}"/2_filtered.fasta
# The 32bit version of usearch will not accept an input file greater than 4GB. The 64bit usearch is $900. Thus, for now:
  echo $(date +%H:%M) "Decompressing merged reads..."
  "${ZIPPER}" -d "${MERGED_READS}".gz

  echo "Calculating merged file size..."
  INFILE_SIZE=$(stat "${MERGED_READS}" | awk '{ print $8 }')
  if [ ${INFILE_SIZE} -gt 4000000000 ]; then
  # Must first check the number of reads. If odd, file must be split so as not to split the middle read's sequence from its quality score.
    echo $(date +%H:%M) "Splitting large input file for quality filtering..."
    LINES_MERGED=$(wc -l < "${MERGED_READS}")
    READS_MERGED=$(( LINES_MERGED / 4 ))
    HALF_LINES=$((LINES_MERGED / 2))
    if [ $((READS_MERGED%2)) -eq 0 ]; then
      head -n ${HALF_LINES} "${MERGED_READS}" > "${MERGED_READS%.*}"_A.fastq
      tail -n ${HALF_LINES} "${MERGED_READS}" > "${MERGED_READS%.*}"_B.fastq
    else
      head -n $(( HALF_LINES + 2 )) "${MERGED_READS}" > "${MERGED_READS%.*}"_A.fastq
      tail -n $(( HALF_LINES - 2 )) "${MERGED_READS}" > "${MERGED_READS%.*}"_B.fastq
    fi
    echo  $(date +%H:%M) "usearch is performing quality control on merged reads..."
    usearch -fastq_filter "${MERGED_READS%.*}"_A.fastq -fastq_maxee "${Max_Expected_Errors}" -fastaout "${LIB_OUTPUT_DIR}"/2_filtered_A.fasta
    usearch -fastq_filter "${MERGED_READS%.*}"_B.fastq -fastq_maxee "${Max_Expected_Errors}" -fastaout "${LIB_OUTPUT_DIR}"/2_filtered_B.fasta
    cat "${LIB_OUTPUT_DIR}"/2_filtered_A.fasta "${LIB_OUTPUT_DIR}"/2_filtered_B.fasta > "${FILTERED_OUTPUT%.*}"_linebreaks.fasta
  else
    echo  $(date +%H:%M) "usearch is performing quality control on merged reads..."
    usearch -fastq_filter "${MERGED_READS}" -fastq_maxee "${Max_Expected_Errors}" -fastaout "${FILTERED_OUTPUT%.*}"_linebreaks.fasta
  fi

  # Remove annoying usearch linebreaks at 80 characters
  echo  $(date +%H:%M) "removing fasta linebreaks..."
  awk '/^>/{print (NR==1)?$0:"\n"$0;next}{printf "%s", $0}END{print ""}' "${FILTERED_OUTPUT%.*}"_linebreaks.fasta > "${FILTERED_OUTPUT}"
  # remove file with linebreaks
  rm "${FILTERED_OUTPUT%.*}"_linebreaks.fasta
fi
