#!/bin/bash

PRIMER_TAGS='/Users/threeprime/Documents/Data/IlluminaData/16S/tags_16S.txt'
TAGS=$(tr '\n' ' ' < "${PRIMER_TAGS}" )
ANALYSIS_DIR='/Users/threeprime/Documents/Data/IlluminaData/16S/20141020/demultest'
DEMULTIPLEX_INPUT='/Users/threeprime/Documents/Data/IlluminaData/16S/20141020/demultest/2_filtered_renamed.fasta'

mkdir "${ANALYSIS_DIR}"/demultiplexed

echo "Demultiplexing: finding left tag (started at $(date +%H:%M))"
for TAG_SEQ in $TAGS; do
(	TAG_DIR="${ANALYSIS_DIR}"/demultiplexed/tag_"${TAG_SEQ}"
	mkdir "${TAG_DIR}"
	awk 'gsub(/^.{0,9}'"$TAG_SEQ"'/,"") {if (a && a !~ /^.{0,9}'"$TAG_SEQ"'/) print a; print } {a = $0}' "${DEMULTIPLEX_INPUT}" > "${TAG_DIR}"/1_tagL_removed.fasta ) &
# 	awk 'gsub(/^.{0,9}'"$TAG_SEQ"'/,"") {if (a && a !~ /^.{0,9}'"$TAG_SEQ"'/) print a; print } {a = $0}' "${DEMULTIPLEX_INPUT}" > "${TAG_DIR}"/1_tagL_present.fasta ) &
# 	awk '/^.{0,9}'"$TAG_SEQ"'/ {line = $0} {print}' "${DEMULTIPLEX_INPUT}" > "${TAG_DIR}"/1_tagL_present.fasta ) &
done

wait

# Remove tags from left side of read
# echo "Demultiplexing: removing left tag (started at $(date +%H:%M))"
# for TAG_SEQ in $TAGS; do
# (	TAG_DIR="${ANALYSIS_DIR}"/demultiplexed/tag_"${TAG_SEQ}"
# 	sed -E 's/^.{0,9}'"${TAG_SEQ}"'//' "${TAG_DIR}"/1_tagL_present.fasta > "${TAG_DIR}"/2_tagL_removed.fasta ) &
# done
# 
# wait

echo "Demultiplexing: trying to find and remove left tag (started at $(date +%H:%M))"
for TAG_SEQ in $TAGS; do
(	TAG_DIR="${ANALYSIS_DIR}"/demultiplexed/tag_"${TAG_SEQ}"
	TAG_RC=$( echo ${TAG_SEQ} | tr "[ATGCatgcNn]" "[TACGtacgNn]" | rev )
	awk 'gsub(/'"$TAG_RC"'.{0,9}$/,"") {if (a && a !~ /'"$TAG_RC"'.{0,9}$/) print a "tag_""'"$TAG_SEQ"'"; print } {a = $0}' "${TAG_DIR}"/1_tagL_removed.fasta > "${TAG_DIR}"/2_oneswoop.fasta ) &
# 	awk 'gsub(/^.{0,9}'"$TAG_SEQ"'/,"") {if (a && a !~ /^.{0,9}'"$TAG_SEQ"'/) print a; print } {a = $0}' "${DEMULTIPLEX_INPUT}" > "${TAG_DIR}"/1_tagL_present.fasta ) &
# 	awk '/^.{0,9}'"$TAG_SEQ"'/ {line = $0} {print}' "${DEMULTIPLEX_INPUT}" > "${TAG_DIR}"/1_tagL_present.fasta ) &
done

wait

# Identify reads containing tags towards the right side of the read
# echo "Demultiplexing: finding right tag (started at $(date +%H:%M))"
# for TAG_SEQ in $TAGS; do
# (	TAG_DIR="${ANALYSIS_DIR}"/demultiplexed/tag_"${TAG_SEQ}"
# 	TAG_RC=$( echo ${TAG_SEQ} | tr "[ATGCatgcNn]" "[TACGtacgNn]" | rev )
# 	grep -E "${TAG_RC}.{0,9}$" -B 1 "${TAG_DIR}"/1_tagL_removed.fasta | grep -v -- "^--$"  > "${TAG_DIR}"/3_tagR_present.fasta ) &
# done
# 
# wait
# 
# echo "Demultiplexing: removing right tag (started at $(date +%H:%M))"
# for TAG_SEQ in $TAGS; do
# (	TAG_DIR="${ANALYSIS_DIR}"/demultiplexed/tag_"${TAG_SEQ}"
# 	TAG_RC=$( echo ${TAG_SEQ} | tr "[ATGCatgcNn]" "[TACGtacgNn]" | rev )
# 	sed -E 's/'"${TAG_RC}"'.{0,9}$//' "${TAG_DIR}"/3_tagR_present.fasta > "${TAG_DIR}"/4_tagR_removed.fasta ) &
# done
# 
# wait
# "${TAG_DIR}"/2_oneswoop.fasta
# echo "Demultiplexing: adding tag sequence to sequenceID (started at $(date +%H:%M))"
# for TAG_SEQ in $TAGS; do
# (	TAG_DIR="${ANALYSIS_DIR}"/demultiplexed/tag_"${TAG_SEQ}"
# 	awk '/^>/ {print $0 "tag_""'"$TAG_SEQ"'"}' "${TAG_DIR}"/2_oneswoop.fasta > "${TAG_DIR}"/5_demult.fasta) &
# # 	awk '/^.{0,9}'"$TAG_SEQ"'/{if (a && a !~ /^.{0,9}'"$TAG_SEQ"'/) print a "tag_""'"$TAG_SEQ"'"; print} {a=$0}' "${DEMULTIPLEX_INPUT}" > "${TAG_DIR}"/1_tagL_present.fasta ) &
# done
# 
# wait
# 
