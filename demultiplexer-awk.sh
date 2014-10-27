#!/bin/bash

cd '/Users/threeprime/Documents/Data/IlluminaData/12S/20140930/Analysis_20141023_0908'

PRIMER_TAGS='/Users/threeprime/Documents/Data/IlluminaData/12S/12s_Tags.txt'

mkdir '/Users/threeprime/Documents/Data/IlluminaData/12S/20140930/Analysis_20141023_0908/demultiplex-awk'

N_TAGS=$( wc -l < "${PRIMER_TAGS}" )

for (( i=1; i<=${N_TAGS}; i++ ));

do
  TAG=$( sed -n $((i * 2))p '/Users/threeprime/Documents/Data/IlluminaData/12S/20140930/Analysis_20141023_0908/tags.fasta' )

  mkdir /Users/threeprime/Documents/Data/IlluminaData/12S/20140930/Analysis_20141023_0908/demultiplex-awk/tag_"${i}"

  CURRENT_DIR=/Users/threeprime/Documents/Data/IlluminaData/12S/20140930/Analysis_20141023_0908/demultiplex-awk/tag_"${i}"
  awk '/^.{0,9}'"$TAG"'/{if (a && a !~ /^.{0,9}'"$TAG"'/) print a,"TAG:""'"$TAG"'"; print} {a=$0}' 3_no_homopolymers.fasta > "${CURRENT_DIR}"/1_no5primetag.fasta

done
