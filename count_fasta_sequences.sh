#!/bin/bash

# This script will find any file named TARGET_FILE in any subdirectory of directory MY_DIR, and count the number of lines containing '>' (i.e. the number of sequences in a FASTA file).

MY_DIR='/Users/threeprime/Documents/GoogleDrive/Illumina_Data/12S/run_20140930/Analysis_20141028_1812_JP/demultiplexed'
TARGET_FILE='3prime_tag_rm.fasta'
OUTPUT_FILE='/Users/threeprime/Desktop/output.txt'

FILES=$( find "${MY_DIR}" -type f -name "${TARGET_FILE}" )

# grep option -c counts selected lines; -H prints filename with output
for FILE in $FILES; do
	grep -c '>' "${FILE}" >> "${OUTPUT_FILE}"
done


# Total sequenced reads
wc -l "${READ1}"
# wc option < before file path suppresses printing of file name/path
LINES_READ1=$( wc -l "${READ1}" )


ANALYSIS_DIR

# Merging of paired-end reads ("successfully merged ?(paired-end) reads"
MERGED_READS

# Fastq quality filtering
FILTERED_OUTPUT

LINES_MERGED=$(wc -l < "${MERGED_READS}")
READS_MERGED=$(( LINES_MERGED / 4 ))

echo $(( 100 / 3 ))
SUMMARY_TABLE=""


echo "collapse rank='$COLLAPSE_RANK1';" >> "${SUMMARY_TABLE}"


# homopolymers
"${ANALYSIS_DIR}"/homopolymer_line_numbers.txt

#
#
# Primer and tag identification
# LEFT TAG
for TAG_SEQ in $TAGS; do
	(	TAG_DIR="${ANALYSIS_DIR}"/demultiplexed/tag_"${TAG_SEQ}"
	mkdir "${TAG_DIR}"
	awk 'gsub(/^.{0,9}'"$TAG_SEQ"'/,"") {if (a && a !~ /^.{0,9}'"$TAG_SEQ"'/) print a; print} {a=$0}' "${DEMULTIPLEX_INPUT}" > "${TAG_DIR}"/1_tagL_removed.fasta ) &
	# awk '/^.{0,9}'"$TAG_SEQ"'/{if (a && a !~ /^.{0,9}'"$TAG_SEQ"'/) print a; print} {a=$0}' "${DEMULTIPLEX_INPUT}" > "${TAG_DIR}"/1_tagL_present.fasta ) &
done

wait

# RIGHT TAG
for TAG_SEQ in $TAGS; do
	(	TAG_DIR="${ANALYSIS_DIR}"/demultiplexed/tag_"${TAG_SEQ}"
	TAG_RC=$( echo ${TAG_SEQ} | tr "[ATGCatgcNn]" "[TACGtacgNn]" | rev )
	awk 'gsub(/'"$TAG_RC"'.{0,9}$/,"") {if (a && a !~ /'"$TAG_RC"'.{0,9}$/) print a "tag_""'"$TAG_SEQ"'"; print } {a = $0}' "${TAG_DIR}"/1_tagL_removed.fasta > "${TAG_DIR}"/2_notags.fasta ) &
done

wait

# PRIMERS
echo "Removing primers..."
for TAG_SEQ in $TAGS; do
	TAG_DIR="${ANALYSIS_DIR}"/demultiplexed/tag_"${TAG_SEQ}"
	# REMOVE PRIMER SEQUENCES
	# Remove PRIMER1 from the beginning of the reads. NOTE cutadapt1.7+ will accept ambiguities in primers.
	cutadapt -g ^"${PRIMER1}" -e "${PRIMER_MISMATCH_PROPORTION}" -m "${LENGTH_ROI_HALF}" --discard-untrimmed "${TAG_DIR}"/2_notags.fasta > "${TAG_DIR}"/5_primerL1_removed.fasta
	cutadapt -g ^"${PRIMER2}" -e "${PRIMER_MISMATCH_PROPORTION}" -m "${LENGTH_ROI_HALF}" --discard-untrimmed "${TAG_DIR}"/2_notags.fasta > "${TAG_DIR}"/5_primerL2_removed.fasta
	# Remove the primer on the other end of the reads by reverse-complementing the files and then trimming PRIMER1 and PRIMER2 from the left side.
	# NOTE cutadapt1.7 will account for anchoring these to the end of the read with $
	seqtk seq -r "${TAG_DIR}"/5_primerL1_removed.fasta | cutadapt -g ^"${PRIMER2}" -e "${PRIMER_MISMATCH_PROPORTION}" -m "${LENGTH_ROI_HALF}" --discard-untrimmed - > "${TAG_DIR}"/6_primerR1_removed.fasta
	seqtk seq -r "${TAG_DIR}"/5_primerL2_removed.fasta | cutadapt -g ^"${PRIMER1}" -e "${PRIMER_MISMATCH_PROPORTION}" -m "${LENGTH_ROI_HALF}" --discard-untrimmed - > "${TAG_DIR}"/6_primerR2_removed.fasta
	seqtk seq -r "${TAG_DIR}"/6_primerR1_removed.fasta > "${TAG_DIR}"/6_primerR1_removedRC.fasta
	cat "${TAG_DIR}"/6_primerR1_removedRC.fasta "${TAG_DIR}"/6_primerR2_removed.fasta > "${TAG_DIR}"/7_no_primers.fasta
done



# Dereplication 
if [ "$CONCATENATE_SAMPLES" = "YES" ]; then
	DEREP_INPUT="${ANALYSIS_DIR}"/concatenated/1_demult_concat.fasta



# Singleton removal
#
# BLAST (percent identity ≥ 98% to 12s mtDNA database)
# Removal of low-frequency noise
# Removal of contaminant reads present in negative controls
# Removal of taxa present in less than 2 out of 3 replicates
# Assignment of false-positive taxa and/or reclassification at lower taxonomic rank
