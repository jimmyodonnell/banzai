#!/bin/bash


# Compute read length (read length is calculated from length of the second line in READ1 file)
# LENGTH_READ=$( sed '2q;d' "${READ1}" | awk '{ print length }' )


# Get the directory containing the READ1 file and assign it to variable READ_DIR.
# READ_DIR="${READ1%/*}"


# PRIMER1RC=$( seqtk seq -r "${PRIMER_FILE}" | awk 'NR==2' )
# PRIMER2RC=$( seqtk seq -r "${PRIMER_FILE}" | awk 'NR==4' )

# DEMULTIPLEXING

# N_TAGS=$( wc -l < "${PRIMER_TAGS}" )

# Write a file of sequence names to make a tag fasta file (necessary for reverse complementing)
# for i in `seq ${N_TAGS}`; do echo \>tag"$i"; done > "${CURRENT_LIB}"/tag_names.txt
# Alternately paste those names and the sequences to make a tag fasta file.
# paste -d"\n" "${CURRENT_LIB}"/tag_names.txt "${PRIMER_TAGS}" > "${CURRENT_LIB}"/tags.fasta
# Reverse complement the tags
# seqtk seq -r "${CURRENT_LIB}"/tags.fasta > "${CURRENT_LIB}"/tags_RC.fasta


# The old loop: Start the loop, do one loop for each of the number of lines in the tag file.
# for (( i=1; i<=${N_TAGS}; i++ ));

# Remove tags from left side of read
# echo "Demultiplexing: removing left tag (started at $(date +%H:%M))"
# for TAG_SEQ in $TAGS; do
# (	TAG_DIR="${CURRENT_LIB}"/demultiplexed/tag_"${TAG_SEQ}"
# 	sed -E 's/^.{0,9}'"${TAG_SEQ}"'//' "${TAG_DIR}"/1_tagL_present.fasta > "${TAG_DIR}"/2_tagL_removed.fasta ) &
# done
#
# wait
#
# Identify reads containing tags towards the right side of the read
# echo "Demultiplexing: finding right tag (started at $(date +%H:%M))"
# for TAG_SEQ in $TAGS; do
# (	TAG_DIR="${CURRENT_LIB}"/demultiplexed/tag_"${TAG_SEQ}"
# 	TAG_RC=$( echo ${TAG_SEQ} | tr "[ATGCatgcNn]" "[TACGtacgNn]" | rev )
# 	grep -E "${TAG_RC}.{0,9}$" -B 1 "${TAG_DIR}"/1_tagL_removed.fasta | grep -v -- "^--$"  > "${TAG_DIR}"/3_tagR_present.fasta ) &
# done
#
# wait
#

# echo "Demultiplexing: adding tag sequence to sequenceID (started at $(date +%H:%M))"
# for TAG_SEQ in $TAGS; do
# (	TAG_DIR="${CURRENT_LIB}"/demultiplexed/tag_"${TAG_SEQ}"
# 	awk '/^.{0,9}'"$TAG_SEQ"'/{if (a && a !~ /^.{0,9}'"$TAG_SEQ"'/) print a "tag_""'"$TAG_SEQ"'"; print} {a=$0}' "${DEMULTIPLEX_INPUT}" > "${TAG_DIR}"/1_tagL_present.fasta ) &
# done

# Assign the current tag to to variable TAG, and its reverse complement to TAG_RC
# TAG=$( sed -n $((i * 2))p "${CURRENT_LIB}"/tags.fasta )
# TAG_RC=$( sed -n $((i * 2))p "${CURRENT_LIB}"/tags_RC.fasta )

# Create a directory for the tag
# mkdir "${CURRENT_LIB}"/demultiplexed/tag_"${i}"

# Make a variable (CURRENT_DIR) with the current tag's directory for ease of reading and writing.
# CURRENT_DIR="${CURRENT_LIB}"/demultiplexed/tag_"${i}"

# REMOVE TAG SEQUENCES
# remove the tag from the beginning of the sequence (5' end) in it's current orientation
# cutadapt -g ^NNN"${TAG}" -e 0 --discard-untrimmed "${DEMULTIPLEX_INPUT}" > "${CURRENT_DIR}"/5prime_tag_rm.fasta

# Need to first grep lines containing pattern first, THEN following sed command with remove them
# grep -E "${TAG_RC}.{0,9}$" -B 1 "${CURRENT_DIR}"/5prime_tag_rm.fasta | grep -v -- "^--$"  > "${CURRENT_DIR}"/3_prime_tagged.fasta
# Turn the following line on to write chimaeras to a file.
# grep -E -v "${TAG_RC}.{0,9}$" "${CURRENT_DIR}"/5prime_tag_rm.fasta > "${CURRENT_DIR}"/chimaeras.fasta
# This sed command looks really f***ing ugly; but I'm pretty sure it works.
# sed -E 's/'"${TAG_RC}"'.{0,9}$//' "${CURRENT_DIR}"/3_prime_tagged.fasta > "${CURRENT_DIR}"/3prime_tag_rm.fasta


# echo "$SCRIPT_DIR/analyses_prelim.R"
# echo "${OUTPUT_PDF}"
# echo "${DEREP_INPUT%/*}"/dups.csv
# echo "${SEQUENCING_POOL_DATA}"
