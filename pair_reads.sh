#!/usr/bin/env bash



MERGED_READS_PREFIX="${LIB_OUTPUT_DIR}"/1_merged
MERGED_READS="${LIB_OUTPUT_DIR}"/1_merged.assembled.fastq

pear \
    -f "${READ1}" \
    -r "${READ2}" \
    -o "${MERGED_READS_PREFIX}" \
    -v $MINOVERLAP \
    -m $ASSMAX \
    -n $ASSMIN \
    -t $TRIMMIN \
    -q $Quality_Threshold \
    -u $UNCALLEDMAX \
    -g $TEST \
    -p $PVALUE \
    -s $SCORING \
    -j $n_threads
