#!/bin/bash

PRIMER_FILE='/Users/threeprime/Documents/Data/IlluminaData/16S/primers_16S.fasta'
PRIMER_MISMATCH_PROPORTION="0.10"

LENGTH_ROI_HALF="69"

PRIMER1=$( awk 'NR==2' "${PRIMER_FILE}" )
PRIMER2=$( awk 'NR==4' "${PRIMER_FILE}" )
# PRIMER1RC=$( seqtk seq -r "${PRIMER_FILE}" | awk 'NR==2' )
# PRIMER2RC=$( seqtk seq -r "${PRIMER_FILE}" | awk 'NR==4' )
PRIMER1RC=$( echo ${PRIMER1} | tr "[ABCDGHMNRSTUVWXYabcdghmnrstuvwxy]" "[TVGHCDKNYSAABWXRtvghcdknysaabwxr]" | rev )
PRIMER2RC=$( echo ${PRIMER2} | tr "[ABCDGHMNRSTUVWXYabcdghmnrstuvwxy]" "[TVGHCDKNYSAABWXRtvghcdknysaabwxr]" | rev )

# REMOVE THE FOLLOWING ONCE WE ALL HAVE CUTADAPT >1.7 UP AND RUNNING
# Take ambiguities out of primers. The sed command says "turn any character that's not A, T, C, or G, and replace it with N.
PRIMER1_NON=$( echo $PRIMER1 | sed "s/[^ATCG]/N/g" )
PRIMER2_NON=$( echo $PRIMER2 | sed "s/[^ATCG]/N/g" )
PRIMER1RC_NON=$( echo $PRIMER1RC | sed "s/[^ATCG]/N/g" )
PRIMER2RC_NON=$( echo $PRIMER2RC | sed "s/[^ATCG]/N/g" )

DEMULTIPLEXED_DIR='/Users/threeprime/Documents/Data/IlluminaData/16S/20141020/Analysis_20141110_1117/demultiplexed'
# for TAG_DIR in $( ls "${DEMULTIPLEXED_DIR}" ); do
# 	seqtk seq -r ${DEMULTIPLEXED_DIR}/${TAG_DIR}/5_primerL1_removed.fasta > ${DEMULTIPLEXED_DIR}/${TAG_DIR}/5_primerL1_removedRC.fasta
# 	seqtk seq -r ${DEMULTIPLEXED_DIR}/${TAG_DIR}/5_primerL2_removed.fasta > ${DEMULTIPLEXED_DIR}/${TAG_DIR}/5_primerL2_removedRC.fasta	
# done

for TAG_DIR in $( ls "${DEMULTIPLEXED_DIR}" ); do
	cutadapt -g ^"${PRIMER2_NON}" -e "${PRIMER_MISMATCH_PROPORTION}" -m "${LENGTH_ROI_HALF}" --discard-untrimmed ${DEMULTIPLEXED_DIR}/${TAG_DIR}/5_primerL1_removedRC.fasta > ${DEMULTIPLEXED_DIR}/"${TAG_DIR}"/5_primerboth_removed1RC.fasta
	cutadapt -g ^"${PRIMER1_NON}" -e "${PRIMER_MISMATCH_PROPORTION}" -m "${LENGTH_ROI_HALF}" --discard-untrimmed ${DEMULTIPLEXED_DIR}/${TAG_DIR}/5_primerL2_removedRC.fasta > ${DEMULTIPLEXED_DIR}/"${TAG_DIR}"/5_primerboth_removed2RC.fasta
	seqtk seq -r ${DEMULTIPLEXED_DIR}/"${TAG_DIR}"/5_primerboth_removed1RC.fasta > ${DEMULTIPLEXED_DIR}/"${TAG_DIR}"/5_primerboth_removed1.fasta
	cat ${DEMULTIPLEXED_DIR}/"${TAG_DIR}"/5_primerboth_removed1.fasta ${DEMULTIPLEXED_DIR}/"${TAG_DIR}"/5_primerboth_removed2RC.fasta > ${DEMULTIPLEXED_DIR}/"${TAG_DIR}"/7_noprimer_same_orient.fasta
done
