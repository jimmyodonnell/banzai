#!/bin/bash

PRIMER_TAGS='/Users/threeprime/Documents/Data/IlluminaData/16S/tags_16S.txt'
TAGS=$(tr '\n' ' ' < "${PRIMER_TAGS}" )
ANALYSIS_DIR='/Users/threeprime/Documents/Data/IlluminaData/16S/20141020/demultest'
PRIMER_FILE='/Users/threeprime/Documents/Data/IlluminaData/16S/primers_16S.fasta'
PRIMER_MISMATCH_PROPORTION="0.10"
LENGTH_ROI_HALF="70"
PRIMER1=$( awk 'NR==2' "${PRIMER_FILE}" )
PRIMER2=$( awk 'NR==4' "${PRIMER_FILE}" )
PRIMER1_NON=$( echo $PRIMER1 | sed "s/[^ATCG]/N/g" )
PRIMER2_NON=$( echo $PRIMER2 | sed "s/[^ATCG]/N/g" )

echo "Removing primers..."
for TAG_SEQ in $TAGS; do
(	TAG_DIR="${ANALYSIS_DIR}"/oneswoop2/tag_"${TAG_SEQ}"
	# REMOVE PRIMER SEQUENCES
	# Remove PRIMER1 from the beginning of the reads. NOTE cutadapt1.7+ will accept ambiguities in primers.
	cutadapt -g ^"${PRIMER1_NON}" -e "${PRIMER_MISMATCH_PROPORTION}" -m "${LENGTH_ROI_HALF}" --discard-untrimmed "${TAG_DIR}"/2_oneswoop.fasta > "${TAG_DIR}"/5_primerL1_removed.fasta
	cutadapt -g ^"${PRIMER2_NON}" -e "${PRIMER_MISMATCH_PROPORTION}" -m "${LENGTH_ROI_HALF}" --discard-untrimmed "${TAG_DIR}"/2_oneswoop.fasta > "${TAG_DIR}"/5_primerL2_removed.fasta
	# Remove the primer on the other end of the reads by reverse-complementing the files and then trimming PRIMER1 and PRIMER2 from the left side.
	# NOTE cutadapt1.7 will account for anchoring these to the end of the read with $
	seqtk seq -r "${TAG_DIR}"/5_primerL1_removed.fasta | cutadapt -g ^"${PRIMER2_NON}" -e "${PRIMER_MISMATCH_PROPORTION}" -m "${LENGTH_ROI_HALF}" --discard-untrimmed - > "${TAG_DIR}"/6_primerR1_removed.fasta
	seqtk seq -r "${TAG_DIR}"/5_primerL2_removed.fasta | cutadapt -g ^"${PRIMER1_NON}" -e "${PRIMER_MISMATCH_PROPORTION}" -m "${LENGTH_ROI_HALF}" --discard-untrimmed - > "${TAG_DIR}"/6_primerR2_removed.fasta
	seqtk seq -r "${TAG_DIR}"/6_primerR1_removed.fasta > "${TAG_DIR}"/6_primerR1_removedRC.fasta
	cat "${TAG_DIR}"/6_primerR1_removedRC.fasta "${TAG_DIR}"/6_primerR2_removed.fasta > "${TAG_DIR}"/7_no_primers.fasta ) &
done

wait
