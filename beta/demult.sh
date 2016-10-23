#!/usr/bin/env bash

################################################################################
# DEMULTIPLEXING (awk)
################################################################################
# input is a single fasta file.
DEMULTIPLEX_INPUT="${FILTERED_OUTPUT}"

# make a directory to put all the demultiplexed files in
DEMULTIPLEXED_DIR="${DEMULTIPLEX_INPUT%/*}"/demultiplexed
mkdir "${DEMULTIPLEXED_DIR}"

# Copy sequences to fasta files into separate directories based on tag sequence on left side of read
# TODO test for speed against removing the tag while finding it: wrap first tag regex in gsub(/pattern/,""):  awk 'gsub(/^.{0,9}'"$TAG_SEQ"'/,""){if . . .
# 20150522 changed {0,9} to {3} to eliminate flexibility (that could result in a read being assigned to >1 sample)
# awk '/^.{0,9}'"$TAG_SEQ"'/{if (a && a !~ /^.{0,9}'"$TAG_SEQ"'/) print a; print} {a=$0}' "${DEMULTIPLEX_INPUT}" > "${TAG_DIR}"/1_tagL_present.fasta ) &


TAG_Ns=3

if [ "${TAG_Ns}" = 0 ]; then
    demult_comm_L='gsub(/^'"$TAG_SEQ"'/,"stuff")
                        {if (a && a !~ /^'"$TAG_SEQ"'/)
                        print a;
                        print
                        } {a=$0}'
    demult_comm_R='gsub(/'"$TAG_RC"'$/,"")
                        {if (a && a !~ /'"$TAG_RC"'$/)
                        print a "tag_""'"$TAG_SEQ"'";
                        print
                        } {a = $0}'
else
    demult_comm_L='gsub(/^.{'"${TAG_Ns}"'}'"$TAG_SEQ"'/,"stuff") {
                        if (a && a !~ /^.{'"${TAG_Ns}"'}'"$TAG_SEQ"'/)
                        print a;
                        print
                        } {a=$0}'
    demult_comm_R='gsub(/'"$TAG_RC"'.{'"${TAG_Ns}"'}$/,"") {
                        if (a && a !~ /'"$TAG_RC"'.{'"${TAG_Ns}"'}$/)
                        print a "tag_""'"$TAG_SEQ"'";
                        print
                        } {a = $0}'
fi

echo $(date +%H:%M) "Demultiplexing: removing tags and adding to sequence ID in library" "${CURRENT_LIB##*/}""..."
for TAG_SEQ in $TAGS; do
    (	TAG_DIR="${DEMULTIPLEXED_DIR}"/tag_"${TAG_SEQ}"
    mkdir "${TAG_DIR}"
    demult_file_L="${TAG_DIR}"/1_tagL_removed.fasta
    demult_file_R="${TAG_DIR}"/2_notags.fasta

    # Left side tag
    awk "${demult_comm_L}" "${DEMULTIPLEX_INPUT}" > "${demult_file_L}"

    # Right side tag
    TAG_RC=$( echo ${TAG_SEQ} | tr "[ATGCatgcNn]" "[TACGtacgNn]" | rev )
    awk  "${demult_comm_R}" "${demult_file_L}" > "${demult_file_R}"

    echo "${CURRENT_LIB##*/}" "${TAG_SEQ}" $(wc -l "${demult_file_L}" | \
        awk '{ print ($1/2) }') $(wc -l "${demult_file_R}" | \
        awk '{ print ($1/2)}') >> "${TAG_COUNT}" ) &

done

wait
