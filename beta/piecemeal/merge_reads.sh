#!/usr/bin/env bash

##############################################################################
# MERGE PAIRED-END READS AND QUALITY FILTER (PEAR)
##############################################################################

LENGTH_READ=$( head -n 100000 "${READ1}" | awk '{print length($0);}' | sort -nr | uniq | head -n 1 )

if [ "${calculate_merge_length}" = "YES" ]; then
  ##############################################################################
  # CALCULATE EXPECTED AND MINIMUM OVERLAP OF PAIRED END SEQUENCES
  ##############################################################################
  OVERLAP_EXPECTED=$(($LENGTH_FRAG - (2 * ($LENGTH_FRAG - $LENGTH_READ) ) ))
  MINOVERLAP=$(( $OVERLAP_EXPECTED / 2 ))
  ##############################################################################
  # CALCULATE MAXIMUM AND MINIMUM LENGTH OF MERGED READS
  ##############################################################################
  ASSMAX=$(( $LENGTH_FRAG + 50 ))
  ASSMIN=$(( $LENGTH_FRAG - 50 ))
else
  MINOVERLAP="${minimum_overlap}"
  ASSMAX="${assembled_max}"
  ASSMIN="${assembled_min}"
fi

if [ "$ALREADY_PEARED" = "YES" ]; then
  MERGED_READS="$PEAR_OUTPUT"
  echo "Paired reads have already been merged."
else
  echo $(date +%H:%M) "Merging reads in library" "${CURRENT_LIB##*/}""..."
  MERGED_READS_PREFIX="${LIB_OUTPUT_DIR}"/1_merged
  MERGED_READS="${LIB_OUTPUT_DIR}"/1_merged.assembled.fastq
  pear \
    --forward-fastq "${READ1}" \
    --reverse-fastq "${READ2}" \
    --output "${MERGED_READS_PREFIX}" \
    -v $MINOVERLAP \
    -m $ASSMAX \
    -n $ASSMIN \
    -t $min_seq_length \
    -q $Quality_Threshold \
    -u $UNCALLEDMAX \
    -g $TEST \
    -p $PVALUE \
    -s $SCORING \
    -j $n_cores

  # check pear output:
  if [[ ! -s "${MERGED_READS}" ]] ; then
      echo 'ERROR: No reads were merged.'
      echo 'Aborting analysis of this library, but will move on to next one.'
      continue
  fi



fi
