#!/bin/bash

# script to perform analysis of primer tagging

ANALYSIS_DIR='/Users/threeprime/Documents/Data/IlluminaData/12S/20140930/Analysis_20141023_0908/demultiplex-awk'

PRIMER_TAGS='/Users/threeprime/Documents/Data/IlluminaData/12S/12s_Tags_error.txt'

N_TAGS=$( wc -l < "${PRIMER_TAGS}" )

for (( i=1; i<=${N_TAGS}; i++ ));

do

  TAG=$( sed -n ${i}p ${PRIMER_TAGS} )

  for N in $(seq 0 1 9);
  do
    VAR+=$( grep -E "^.{$N}${TAG}" "${ANALYSIS_DIR}"/tag_"$i"/1_no5primetag.fasta | wc -l )
    echo "$VAR" >> outputfile
    unset VAR
  done

done





    # echo "$Age,$Gender,$Name,$Street,$Occupation" >> outputfile


# for (( N=0; N<=9; N++ ));
#
# do
#
#   VAR+=$( grep -E "^.{$N}${TAG}" "${ANALYSIS_DIR}"/tag_"$i"/1_no5primetag.fasta | wc -l )
#   echo "Tag found $N bases from start of sequence"
#
# done
