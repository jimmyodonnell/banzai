#!/bin/bash

INFILE='/Users/threeprime/Documents/Data/IlluminaData/12S/20140930/Analysis_20141023_0908/3_no_homopolymers.fasta'

PRIMER_TAGS='/Users/threeprime/Documents/Data/IlluminaData/12S/12s_Tags_error.txt'

NEW_DIR="/Users/threeprime/Documents/Data/IlluminaData/12S/20140930/Analysis_20141023_0908/demultiplex-awk-parallel"

mkdir "${NEW_DIR}"

# N_TAGS=$( wc -l < "${PRIMER_TAGS}" )

TAGS=$(tr '\n' ' ' < "${PRIMER_TAGS}" )

# for (( i=1; i<=${N_TAGS}; i++ ));

for T in $TAGS

do
  # TAG=$( sed -n $((i * 2))p '/Users/threeprime/Documents/Data/IlluminaData/12S/20140930/Analysis_20141023_0908/tags.fasta' )

(  TAG_DIR="${NEW_DIR}"/tag_"${T}"

   mkdir "${TAG_DIR}"

   awk '/^.{0,9}'"$T"'/{if (a && a !~ /^.{0,9}'"$T"'/) print a,"TAG:""'"$T"'"; print} {a=$0}' "${INFILE}" > "${TAG_DIR}"/1_no5primetag.fasta ) &

done
