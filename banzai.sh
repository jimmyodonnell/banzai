#!/usr/bin/env bash

# Pipeline for analysis of MULTIPLEXED Illumina data, a la Jimmy

# TODO An attempt to cause the script to exit if any of the commands returns a non-zero status (i.e. FAILS).
# set -e

################################################################################
# RAW DATA, ANALYSIS PARAMETERS, AND GENERAL SETTINGS
################################################################################

# Define a variable called START_TIME
START_TIME=$(date +%Y%m%d_%H%M)

# This command specifies the path to the directory containing this script
SCRIPT_DIR="$(dirname "$0")"

# Read in the parameter file
source "$SCRIPT_DIR/banzai_params.sh"

# make an analysis directory with starting time timestamp
ANALYSIS_DIR="${ANALYSIS_DIRECTORY}"/Analysis_"${START_TIME}"
mkdir "${ANALYSIS_DIR}"

# Write a log file of output from this script (everything that prints to terminal)
exec > >(tee "${ANALYSIS_DIR}"/logfile.txt) 2>&1
# exec 2>&1

echo "Analysis started at ""${START_TIME}"

# Read in the PEAR parameter files
source "$SCRIPT_DIR/pear_params.sh"

# Detect number of cores on machine; set variable
N_CORES=$(sysctl -n hw.ncpu)
if [ $N_CORES -gt 1 ]; then
	echo "$N_CORES cores detected."
else
	N_CORES=1
	echo "Multiple cores not detected."
fi


# Copy these files into that directory as a verifiable log you can refer back to.
cp "${SCRIPT_DIR}"/banzai.sh "${ANALYSIS_DIR}"/analysis_script.txt
cp "${SCRIPT_DIR}"/banzai_params.sh "${ANALYSIS_DIR}"/analysis_parameters.txt
cp "${SCRIPT_DIR}"/pear_params.sh "${ANALYSIS_DIR}"/pear_parameters.txt

################################################################################
# LOAD MULTIPLEX TAGS
################################################################################
if [ "${READ_TAGS_FROM_SEQUENCING_POOL_DATA}" = "YES" ]; then
	TAG_COL=$(awk -F',' -v TAG_COL_NAME=$TAG_COLUMN_NAME '{for (i=1;i<=NF;i++) if($i ~ TAG_COL_NAME) print i; exit}' $SEQUENCING_POOL_DATA)
	TAGS=$(awk -F',' -v TAGCOL=$TAG_COL 'NR>1 {print $TAGCOL}' $SEQUENCING_POOL_DATA | sort | uniq)
	echo "Multiplex tags read from sequencing pool data."
else
	TAGS=$(tr '\n' ' ' < "${TAG_FILE}" )
	echo "Multiplex tags read from tag file."
fi
# make tag sequences into an array
declare -a TAGS_ARRAY=($TAGS)


################################################################################
# Read in primers and create reverse complements.
################################################################################
if [ "${READ_PRIMERS_FROM_SEQUENCING_POOL_DATA}" = "YES" ]; then
	PRIMER1_COLNUM=$(awk -F',' -v PRIMER1_COL=$PRIMER_1_COLUMN_NAME '{for (i=1;i<=NF;i++) if($i ~ PRIMER1_COL) print i; exit}' $SEQUENCING_POOL_DATA)
	PRIMER2_COLNUM=$(awk -F',' -v PRIMER2_COL=$PRIMER_2_COLUMN_NAME '{for (i=1;i<=NF;i++) if($i ~ PRIMER2_COL) print i; exit}' $SEQUENCING_POOL_DATA)
	PRIMER1=$(awk -F',' -v PRIMER1_COL=$PRIMER1_COLNUM 'NR==2 {print $PRIMER1_COL}' $SEQUENCING_POOL_DATA)
	PRIMER2=$(awk -F',' -v PRIMER2_COL=$PRIMER2_COLNUM 'NR==2 {print $PRIMER2_COL}' $SEQUENCING_POOL_DATA)
	echo "Primers read from sequencing pool data."
else
	PRIMER1=$( awk 'NR==2' "${PRIMER_FILE}" )
	PRIMER2=$( awk 'NR==4' "${PRIMER_FILE}" )
	echo "Primers read from primer file."
fi

PRIMER1RC=$( echo ${PRIMER1} | tr "[ABCDGHMNRSTUVWXYabcdghmnrstuvwxy]" "[TVGHCDKNYSAABWXRtvghcdknysaabwxr]" | rev )
PRIMER2RC=$( echo ${PRIMER2} | tr "[ABCDGHMNRSTUVWXYabcdghmnrstuvwxy]" "[TVGHCDKNYSAABWXRtvghcdknysaabwxr]" | rev )

# Calculate the expected size of the region of interest, given the total size of fragments, and the length of primers and tags
EXTRA_SEQ=${TAGS_ARRAY[0]}${TAGS_ARRAY[0]}$PRIMER1$PRIMER2
LENGTH_ROI=$(( $LENGTH_FRAG - ${#EXTRA_SEQ} ))
LENGTH_ROI_HALF=$(( $LENGTH_ROI / 2 ))

LIBRARY_DIRECTORIES=$( find "$PARENT_DIR" -name '*.fastq' -print0 | xargs -0 -n1 dirname | sort --unique )

# TODO WOULD LIKE TO ADD A LIBRARY NAMES VARIABLE
# LIBRARY_NAMES=$( ls "$PARENT_DIR" )
# declare -a LIBRARY_NAMES_A=($LIBRARY_NAMES)

################################################################################
# BEGIN LOOP TO PERFORM LIBRARY-LEVEL ACTIONS
################################################################################

for CURRENT_LIB in $LIBRARY_DIRECTORIES; do

	# Identify the forward and reverse fastq files.
	READ1=$(find "${CURRENT_LIB}" -name '*R1*fastq')
	READ2=$(find "${CURRENT_LIB}" -name '*R2*fastq')

	LIB_OUTPUT_DIR="${ANALYSIS_DIR}"/${CURRENT_LIB##*/}
	mkdir "${LIB_OUTPUT_DIR}"

	# TODO remove whitespace from sequence labels
	# sed 's/ /_/'

	################################################################################
	# MERGE PAIRED-END READS (PEAR)
	################################################################################
	if [ "$ALREADY_PEARED" = "YES" ]; then
		MERGED_READS="$PEAR_OUTPUT"
		echo "Paired reads have already been merged..."
	else
		MERGED_READS_PREFIX="${LIB_OUTPUT_DIR}"/1_merged
		MERGED_READS="${LIB_OUTPUT_DIR}"/1_merged.assembled.fastq
		pear -f "${READ1}" -r "${READ2}" -o "${MERGED_READS_PREFIX}" -v $MINOVERLAP -m $ASSMAX -n $ASSMIN -t $TRIMMIN -q $QT -u $UNCALLEDMAX -g $TEST -p $PVALUE -s $SCORING -j $THREADS
	fi

	################################################################################
	# QUALITY FILTERING (usearch)
	################################################################################
	# FILTER READS (This is the last step that uses quality scores, so convert to fasta)
	if [ "${ALREADY_FILTERED}" = "YES" ]; then
		echo "Using existing filtered reads in file $FILTERED_OUTPUT"
	else
		FILTERED_OUTPUT="${LIB_OUTPUT_DIR}"/2_filtered.fasta
	# The 32bit version of usearch will not accept an input file greater than 4GB. The 64bit usearch is $900. Thus, for now:
		echo "Calculating merged file size..."
		INFILE_SIZE=$(stat "${MERGED_READS}" | awk '{ print $8 }')
		if [ ${INFILE_SIZE} -gt 4000000000 ]; then
		# Must first check the number of reads. If odd, file must be split so as not to split the middle read's sequence from its quality score.
			echo $(date +%H:%M) "Splitting large input file for quality filtering..."
			LINES_MERGED=$(wc -l < "${MERGED_READS}")
			READS_MERGED=$(( LINES_MERGED / 4 ))
			HALF_LINES=$((LINES_MERGED / 2))
			if [ $((READS_MERGED%2)) -eq 0 ]; then
				head -n ${HALF_LINES} "${MERGED_READS}" > "${MERGED_READS%.*}"_A.fastq
				tail -n ${HALF_LINES} "${MERGED_READS}" > "${MERGED_READS%.*}"_B.fastq
			else
				head -n $(( HALF_LINES + 2 )) "${MERGED_READS}" > "${MERGED_READS%.*}"_A.fastq
				tail -n $(( HALF_LINES - 2 )) "${MERGED_READS}" > "${MERGED_READS%.*}"_B.fastq
			fi
			usearch -fastq_filter "${MERGED_READS%.*}"_A.fastq -fastq_maxee 0.5 -fastq_minlen "${ASSMIN}" -fastaout "${LIB_OUTPUT_DIR}"/2_filtered_A.fasta
			usearch -fastq_filter "${MERGED_READS%.*}"_B.fastq -fastq_maxee 0.5 -fastq_minlen "${ASSMIN}" -fastaout "${LIB_OUTPUT_DIR}"/2_filtered_B.fasta
			cat "${LIB_OUTPUT_DIR}"/2_filtered_A.fasta "${LIB_OUTPUT_DIR}"/2_filtered_B.fasta > "${FILTERED_OUTPUT}"
		else
			echo  $(date +%H:%M) "usearch is performing quality control on merged reads..."
			usearch -fastq_filter "${MERGED_READS}" -fastq_maxee 0.5 -fastq_minlen "${ASSMIN}" -fastaout "${FILTERED_OUTPUT}"_linebreaks
			echo  $(date +%H:%M) "removing fasta linebreaks..."
			awk '/^>/{print (NR==1)?$0:"\n"$0;next}{printf "%s", $0}END{print ""}' "${FILTERED_OUTPUT}"_linebreaks > "${FILTERED_OUTPUT}"
		fi
	fi

	if [ "${RENAME_READS}" = "YES" ]; then
		echo $(date +%H:%M) "Renaming reads..."
		# usearch7:	sed -E "s/ (1|2):N:0:[0-9]/_"${CURRENT_LIB##*/}"_/" "${FILTERED_OUTPUT}" > "${CURRENT_LIB}"/tmp.fasta
		# usearch8, which without warning removes any part of the sequence ID following a space.
		sed -E "s/$/_"${CURRENT_LIB##*/}"_/" "${FILTERED_OUTPUT}" > "${CURRENT_LIB}"/tmp.fasta
		sed -E "s/>([a-zA-Z0-9-]*:){4}/>/" "${CURRENT_LIB}"/tmp.fasta > "${FILTERED_OUTPUT%.*}"_renamed.fasta
		rm "${CURRENT_LIB}"/tmp.fasta
		FILTERED_OUTPUT="${FILTERED_OUTPUT%.*}"_renamed.fasta
	else
		echo "Reads not renamed"
	fi


	################################################################################
	# HOMOPOLYMERS (grep, awk)
	################################################################################
	if [ "${REMOVE_HOMOPOLYMERS}" = "YES" ]; then
		echo $(date +%H:%M) "Removing homopolymers..."
		grep -E -i "(A|T|C|G)\1{$HOMOPOLYMER_MAX,}" "${FILTERED_OUTPUT}" -B 1 -n | cut -f1 -d: | cut -f1 -d- | sed '/^$/d' > "${CURRENT_LIB}"/homopolymer_line_numbers.txt
		if [ -s "${CURRENT_LIB}"/homopolymer_line_numbers.txt ]; then
			awk 'NR==FNR{l[$0];next;} !(FNR in l)' "${CURRENT_LIB}"/homopolymer_line_numbers.txt "${FILTERED_OUTPUT}" > "${CURRENT_LIB}"/3_no_homopolymers.fasta
			awk 'NR==FNR{l[$0];next;} (FNR in l)' "${CURRENT_LIB}"/homopolymer_line_numbers.txt "${FILTERED_OUTPUT}" > "${CURRENT_LIB}"/homopolymeric_reads.fasta
			DEMULTIPLEX_INPUT="${CURRENT_LIB}"/3_no_homopolymers.fasta
		else
			echo "No homopolymers found" > "${CURRENT_LIB}"/3_no_homopolymers.fasta
			DEMULTIPLEX_INPUT="${FILTERED_OUTPUT}"
		fi
	else
		echo "Homopolymers not removed."
		DEMULTIPLEX_INPUT="${FILTERED_OUTPUT}"
	fi

	################################################################################
	# DEMULTIPLEXING (awk)
	################################################################################
	# make a directory to put all the demultiplexed files in
	DEMULTIPLEXED_DIR="${LIB_OUTPUT_DIR}"/demultiplexed
	mkdir "${DEMULTIPLEXED_DIR}"

	# Copy sequences to fasta files into separate directories based on tag sequence on left side of read
	# TODO test for speed against removing the tag while finding it: wrap first tag regex in gsub(/pattern/,""):  awk 'gsub(/^.{0,9}'"$TAG_SEQ"'/,""){if . . .
	echo "Demultiplexing: removing left tag (started at $(date +%H:%M))"
	for TAG_SEQ in $TAGS; do
	(	TAG_DIR="${LIB_OUTPUT_DIR}"/demultiplexed/tag_"${TAG_SEQ}"
		mkdir "${TAG_DIR}"
		awk 'gsub(/^.{0,9}'"$TAG_SEQ"'/,"") {if (a && a !~ /^.{0,9}'"$TAG_SEQ"'/) print a; print} {a=$0}' "${DEMULTIPLEX_INPUT}" > "${TAG_DIR}"/1_tagL_removed.fasta ) &
		# awk '/^.{0,9}'"$TAG_SEQ"'/{if (a && a !~ /^.{0,9}'"$TAG_SEQ"'/) print a; print} {a=$0}' "${DEMULTIPLEX_INPUT}" > "${TAG_DIR}"/1_tagL_present.fasta ) &
	done

	wait

	echo "Demultiplexing: removing right tag and adding tag sequence to sequence ID (started at $(date +%H:%M))"
	for TAG_SEQ in $TAGS; do
	(	TAG_DIR="${LIB_OUTPUT_DIR}"/demultiplexed/tag_"${TAG_SEQ}"
		TAG_RC=$( echo ${TAG_SEQ} | tr "[ATGCatgcNn]" "[TACGtacgNn]" | rev )
		awk 'gsub(/'"$TAG_RC"'.{0,9}$/,"") {if (a && a !~ /'"$TAG_RC"'.{0,9}$/) print a "tag_""'"$TAG_SEQ"'"; print } {a = $0}' "${TAG_DIR}"/1_tagL_removed.fasta > "${TAG_DIR}"/2_notags.fasta ) &
	done

	wait

done
################################################################################
# END LOOP TO PERFORM LIBRARY-LEVEL ACTIONS
################################################################################










# TODO add single if/else for CONCATENATE_SAMPLES: assign directory as appropriate, correct references within loop to be extendable
if [ "$CONCATENATE_SAMPLES" = "YES" ]; then
	# do the stuff here
	WORKING_DIR="${CONCAT_DIR}"
else
	WORKING_DIR="${LIBRARY_DIRECTORIES}"
fi

################################################################################
# CONCATENATE SAMPLES
################################################################################
# TODO could move this first step up above any loops (no else)
if [ "$CONCATENATE_SAMPLES" = "YES" ]; then

	echo $(date +%H:%M) "Concatenating fasta files..."
	CONCAT_DIR="$ANALYSIS_DIR"/all_lib
	mkdir "${CONCAT_DIR}"

	# TODO could move this into above loop after demultiplexing?
	for CURRENT_LIB in $LIBRARY_DIRECTORIES; do

		LIB_OUTPUT_DIR="${ANALYSIS_DIR}"/${CURRENT_LIB##*/}

		for TAG_SEQ in $TAGS; do
			cat "${LIB_OUTPUT_DIR}"/demultiplexed/tag_"${TAG_SEQ}"/2_notags.fasta >> "${CONCAT_DIR}"/1_demult_concat.fasta
		done


	done

	################################################################################
	# PRIMER REMOVAL
	################################################################################
	echo $(date +%H:%M) "Removing primers..."
	# for TAG_SEQ in $TAGS; do
		# TAG_DIR="${ANALYSIS_DIR}"/demultiplexed/tag_"${TAG_SEQ}"
		# Remove PRIMER1 and PRIMER2 from the BEGINNING of the reads. NOTE cutadapt1.7+ will accept ambiguities in primers.
		cutadapt -g ^"${PRIMER1}" -e "${PRIMER_MISMATCH_PROPORTION}" -m "${LENGTH_ROI_HALF}" --discard-untrimmed "${CONCAT_DIR}"/1_demult_concat.fasta > "${CONCAT_DIR}"/5_primerL1_removed.fasta
		cutadapt -g ^"${PRIMER2}" -e "${PRIMER_MISMATCH_PROPORTION}" -m "${LENGTH_ROI_HALF}" --discard-untrimmed "${CONCAT_DIR}"/1_demult_concat.fasta > "${CONCAT_DIR}"/5_primerL2_removed.fasta
		# Remove the primer on the other end of the reads by reverse-complementing the files and then trimming PRIMER1 and PRIMER2 from the left side.
		# NOTE cutadapt1.7 will account for anchoring these to the end of the read with $
		# seqtk seq -r "${CONCAT_DIR}"/5_primerL1_removed.fasta | cutadapt -g ^"${PRIMER2}" -e "${PRIMER_MISMATCH_PROPORTION}" -m "${LENGTH_ROI_HALF}" --discard-untrimmed - > "${CONCAT_DIR}"/6_primerR1_removed.fasta
		cutadapt -a "${PRIMER2RC}"$ -e "${PRIMER_MISMATCH_PROPORTION}" -m "${LENGTH_ROI_HALF}" --discard-untrimmed "${CONCAT_DIR}"/5_primerL1_removed.fasta > "${CONCAT_DIR}"/6_primerR1_removed.fasta
		# seqtk seq -r "${CONCAT_DIR}"/5_primerL2_removed.fasta | cutadapt -g ^"${PRIMER1}" -e "${PRIMER_MISMATCH_PROPORTION}" -m "${LENGTH_ROI_HALF}" --discard-untrimmed - > "${CONCAT_DIR}"/6_primerR2_removed.fasta
		cutadapt -a "${PRIMER1RC}"$ -e "${PRIMER_MISMATCH_PROPORTION}" -m "${LENGTH_ROI_HALF}" --discard-untrimmed "${CONCAT_DIR}"/5_primerL2_removed.fasta > "${CONCAT_DIR}"/6_primerR2_removed.fasta
		seqtk seq -r "${CONCAT_DIR}"/6_primerR1_removed.fasta > "${CONCAT_DIR}"/6_primerR1_removedRC.fasta
		cat "${CONCAT_DIR}"/6_primerR1_removedRC.fasta "${CONCAT_DIR}"/6_primerR2_removed.fasta > "${CONCAT_DIR}"/7_no_primers.fasta
		DEREP_INPUT="${CONCAT_DIR}"/7_no_primers.fasta
	# done

	################################################################################
	# CONSOLIDATE IDENTICAL SEQUENCES
	################################################################################
	echo $(date +%H:%M) "Identifying identical sequences... (python)"
	python "$SCRIPT_DIR/dereplicate_fasta.py" "${DEREP_INPUT}"
	# usearch -derep_fulllength "${DEREP_INPUT}" -sizeout -strand both -uc "${DEREP_INPUT%/*}"/2_derep.uc -output "${DEREP_INPUT%/*}"/2_derep.fasta

	# COUNT DUPLICATES PER READ, REMOVE SINGLETONS
	echo $(date +%H:%M) "Counting duplicates per identical sequence... (awk)"
	awk -F';' '{ if (NF > 2) print NF-1 ";" $0 }' "${DEREP_INPUT}".all | sort -nr | awk -F';' '{ print ">DUP_" NR ";" $0}' > ${DEREP_INPUT%/*}/nosingle.txt

	# COUNT OCCURRENCES PER SAMPLE (LIBRARY + TAG) PER DUPLICATE
	echo $(date +%H:%M) "Consolidating identical sequences per sample... (awk)"
	for CURRENT_LIB in $LIBRARY_DIRECTORIES; do
		for TAG_SEQ in $TAGS; do
			LIB_TAG="${CURRENT_LIB##*/}_tag_${TAG_SEQ}"
			( awk 'BEGIN {print "'$LIB_TAG'" ; FS ="'${LIB_TAG}'" } { print NF -1 }' ${DEREP_INPUT%/*}/nosingle.txt > ${DEREP_INPUT%/*}/"${LIB_TAG}".dup ) &
		done
	done

	wait

	# Write a csv file of the number of occurrences of each duplicate sequence per tag.
	find "${DEREP_INPUT%/*}" -type f -name '*.dup' -exec paste -d, {} \+ | awk '{ print "DUP_" NR-1 "," $0 }' > "${DEREP_INPUT%/*}"/dups.csv
	# delete all of the '.dup' files
	rm ${DEREP_INPUT%/*}/*.dup

	# Write fasta file in order to blast sequences
	echo $(date +%H:%M) "Writing fasta file for duplicate (unBLASTed sequences)"
	awk -F';' '{ print $1 ";size=" $2 ";\n" $3 }' ${DEREP_INPUT%/*}/nosingle.txt > ${DEREP_INPUT%/*}/no_duplicates.fasta

	################################################################################
	# CLUSTER OTUS
	################################################################################
	# Note that identical (duplicate) sequences were consolidated earlier;
	# This step outputs a file (*.uc) that lists, for every sequence, which sequence it clusters with
	if [ "$CLUSTER_OTUS" = "NO" ]; then
		BLAST_INPUT="${DEREP_INPUT%/*}"/no_duplicates.fasta
	else
		echo $(date +%H:%M) "Clustering OTUs..."
		CLUSTER_RADIUS="$(( 100 - ${CLUSTERING_PERCENT} ))"
		UPARSE_OUT="${DEREP_INPUT%/*}"/OTU_uparse.txt
		usearch -cluster_otus "${DEREP_INPUT%/*}"/no_duplicates.fasta -otu_radius_pct "${CLUSTER_RADIUS}" -sizein -sizeout -otus "${DEREP_INPUT%/*}"/9_OTUs_linebreaks.fasta -uparseout "${UPARSE_OUT}"

		# remove the annoying line breaks
		echo $(date +%H:%M) "Removing line breaks within fasta sequences generated by usearch (awk)..."
		awk '/^>/{print (NR==1)?$0:"\n"$0;next}{printf "%s", $0}END{print ""}' "${DEREP_INPUT%/*}"/9_OTUs_linebreaks.fasta > "${DEREP_INPUT%/*}"/9_OTUs.fasta
		# usearch option -notmatched disappeared with version 8
		# awk '/^>/{print (NR==1)?$0:"\n"$0;next}{printf "%s", $0}END{print ""}' "${DEREP_INPUT%/*}"/9_notmatched_linebreaks.fasta > "${DEREP_INPUT%/*}"/9_notmatched.fasta

		rm "${DEREP_INPUT%/*}"/9_OTUs_linebreaks.fasta # "${DEREP_INPUT%/*}"/9_notmatched_linebreaks.fasta

		################################################################################
		# RESOLVE OTUS AND DUPLICATES
		################################################################################
		# Note that this drops sequences determined to be chimeras by usearch
		DUPS_TO_OTUS="${DEREP_INPUT%/*}"/dups_to_otus.csv
		awk -F'[\t;]' 'BEGIN{ print "Query,Match" } { if ($4 == "otu") {print $1 "," $1} else if ($4 == "match") { print $1 "," $7 } else if ($4 == "chimera") { print $1 "," "chimera"} }' "${UPARSE_OUT}" > "${DUPS_TO_OTUS}"

		# Convert duplicate table to OTU table using R script
		Rscript "$SCRIPT_DIR/dup_to_OTU_table.R" "${CONCAT_DIR}"

		BLAST_INPUT="${DEREP_INPUT%/*}"/9_OTUs.fasta
	fi


	################################################################################
	# BLAST CLUSTERS
	################################################################################
	echo $(date +%H:%M) "BLASTing..."
	blastn -query "${BLAST_INPUT}" -db "$BLAST_DB" -num_threads "$N_CORES" -perc_identity "${PERCENT_IDENTITY}" -word_size "${WORD_SIZE}" -evalue "${EVALUE}" -max_target_seqs "${MAXIMUM_MATCHES}" -outfmt 5 -out "${DEREP_INPUT%/*}"/10_BLASTed.xml















else
################################################################################
# DON'T CONCATENATE SAMPLES
################################################################################

	################################################################################
	# PRIMER REMOVAL
	################################################################################
	echo $(date +%H:%M) "Removing primers..."
	for TAG_SEQ in $TAGS; do
		TAG_DIR="${ANALYSIS_DIR}"/demultiplexed/tag_"${TAG_SEQ}"
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

	################################################################################
	# CONSOLIDATE IDENTICAL SEQUENCES
	################################################################################
	for TAG_SEQ in $TAGS; do
		TAG_DIR="${ANALYSIS_DIR}"/demultiplexed/tag_"${TAG_SEQ}"

		DEREP_INPUT="${TAG_DIR}"/7_no_primers.fasta

		# usearch -derep_fulllength "${DEREP_INPUT}" -sizeout -strand both -uc "${TAG_DIR}"/derep.uc -output "${TAG_DIR}"/7_derep.fasta
		echo $(date +%H:%M) "Consolidating identical sequences..."
		python "$SCRIPT_DIR/dereplicate_fasta.py" "${DEREP_INPUT}"

		# REMOVE SINGLETONS
		# usearch -sortbysize "${TAG_DIR}"/7_derep.fasta -minsize 2 -sizein -sizeout -output "${TAG_DIR}"/8_nosingle.fasta
		# COUNT DUPLICATES PER READ, REMOVE SINGLETONS
		awk -F';' '{ if (NF > 2) print NF-1 ";" $0 }' "${DEREP_INPUT}".all | sort -nr | awk -F';' '{ print ">DUP_" NR ";" $0}' > ${DEREP_INPUT%/*}/nosingle.txt

		# count the duplicates
		awk 'BEGIN { FS ="_tag_'${TAG_SEQ}'" } { print NF -1 }' "${DEREP_INPUT%/*}"/nosingle.txt > ${DEREP_INPUT%/*}/"${TAG_SEQ}".dup

		# Write fasta file in order to blast sequences
		awk -F';' '{ print $1 ";size=" $2 ";\n" $3 }' ${DEREP_INPUT%/*}/nosingle.txt > ${DEREP_INPUT%/*}/no_duplicates.fasta

		# CLUSTER SEQUENCES
		if [ "$CLUSTER_OTUS" = "NO" ]; then
			BLAST_INPUT=${DEREP_INPUT%/*}/no_duplicates.fasta
		else
			CLUSTER_RADIUS="$(( 100 - ${CLUSTERING_PERCENT} ))"
			UPARSE_OUT="${DEREP_INPUT%/*}"/OTU_uparse.txt
			usearch -cluster_otus "${DEREP_INPUT%/*}"/nosingle.txt -otu_radius_pct "${CLUSTER_RADIUS}" -sizein -sizeout -otus "${TAG_DIR}"/9_OTUs.fasta -uparseout "${UPARSE_OUT}"
			BLAST_INPUT="${TAG_DIR}"/9_OTUs.fasta
		fi

		# BLAST CLUSTERS
		blastn -query "${BLAST_INPUT}" -db "$BLAST_DB" -num_threads "$N_CORES" -perc_identity "${PERCENT_IDENTITY}" -word_size "${WORD_SIZE}" -evalue "${EVALUE}" -max_target_seqs "${MAXIMUM_MATCHES}" -outfmt 5 -out "${TAG_DIR}"/10_BLASTed.xml
	done
fi



################################################################################
# TAXONOMIC ANNOTATION
################################################################################
if [ "$CONCATENATE_SAMPLES" = "YES" ]; then
	DIRECTORIES="${DEREP_INPUT%/*}"
else
	DIRECTORIES=$( find "${ANALYSIS_DIR}"/demultiplexed -type d -d 1 )
fi

for DIR in "$DIRECTORIES"; do

		# Some POTENTIAL OPTIONS FOR MEGAN EXPORT:
		# {readname_taxonname|readname_taxonid|readname_taxonpath|readname_matches|taxonname_count|taxonpath_count|taxonid_count|taxonname_readname|taxonpath_readname|taxonid_readname}
		# PERFORM COMMON ANCESTOR GROUPING IN MEGAN

		# Specify paths to megan-related files
		BLAST_XML="${DIR}"/10_BLASTed.xml
		MEGAN_COMMAND_FILE="${DIR}"/megan_commands.txt
		MEGAN_RMA_FILE="${DIR}"/meganfile.rma
		MEGAN_SHELL_SCRIPT="${DIR}"/megan_script.sh

		echo "import blastfile='${BLAST_XML}' meganFile='${MEGAN_RMA_FILE}' \
minScore=${MINIMUM_SCORE} \
maxExpected=${MAX_EXPECTED} \
topPercent=${TOP_PERCENT} \
minSupportPercent=${MINIMUM_SUPPORT_PERCENT} \
minSupport=${MINIMUM_SUPPORT} \
minComplexity=${MINIMUM_COMPLEXITY} \
lcapercent=${LCA_PERCENT};" > "${MEGAN_COMMAND_FILE}"
		echo "update;" >> "${MEGAN_COMMAND_FILE}"
		echo "collapse rank='$COLLAPSE_RANK1';" >> "${MEGAN_COMMAND_FILE}"
		echo "update;" >> "${MEGAN_COMMAND_FILE}"
		echo "select nodes=all;" >> "${MEGAN_COMMAND_FILE}"
		echo "export what=DSV format=readname_taxonname separator=comma file=${DIR}/meganout_${COLLAPSE_RANK1}.csv;" >> "${MEGAN_COMMAND_FILE}"
		if [ "$PERFORM_SECONDARY_MEGAN" = "YES" ]; then
			echo "collapse rank='$COLLAPSE_RANK2';" >> "${MEGAN_COMMAND_FILE}"
			echo "update;" >> "${MEGAN_COMMAND_FILE}"
			echo "select nodes=all;" >> "${MEGAN_COMMAND_FILE}"
			echo "export what=DSV format=readname_taxonname separator=comma file=${DIR}/meganout_${COLLAPSE_RANK2}.csv;" >> "${MEGAN_COMMAND_FILE}"
		fi
		echo "quit;" >> "${MEGAN_COMMAND_FILE}"

		echo "#!/bin/bash" > "$MEGAN_SHELL_SCRIPT"
		echo "cd "${megan_exec%/*}"" >> "$MEGAN_SHELL_SCRIPT"
		echo "./"${megan_exec##*/}" -g -E -c ${DIR}/megan_commands.txt" >> "$MEGAN_SHELL_SCRIPT"

		# Run MEGAN
		sh "${DIR}"/megan_script.sh

		# Modify the MEGAN output so that it is a standard CSV file with clusterID, N_reads, and Taxon
		sed 's|;size=|,|' <"${DIR}"/meganout_${COLLAPSE_RANK1}.csv >"${DIR}"/meganout_${COLLAPSE_RANK1}_mod.csv
		sed 's|;size=|,|' <"${DIR}"/meganout_${COLLAPSE_RANK2}.csv >"${DIR}"/meganout_${COLLAPSE_RANK2}_mod.csv

		# Run the R script, passing the current tag directory as the directory to which R will "setwd()"
		Rscript "$SCRIPT_DIR/megan_plotter.R" "${DIR}"

	done


################################################################################
# PRELIMINARY ANALYSES
################################################################################
# Once you have a final CSV file of the number of occurences of each OTU in each sample, run some preliminary analyses in R

OUTPUT_PDF="${ANALYSIS_DIR}"/analysis_results_"${START_TIME}".pdf

echo $(date +%H:%M) "passing args to R..."
Rscript "$SCRIPT_DIR/analyses_prelim.R" "${OUTPUT_PDF}" "${DEREP_INPUT%/*}"/dups.csv "${SEQUENCING_POOL_DATA}"

REMOTE_PDF="${OUTPUT_PDF_DIR}"/analysis_results_"${START_TIME}".pdf
cp "${OUTPUT_PDF}" "${REMOTE_PDF}"


if [ "$PERFORM_CLEANUP" = "YES" ]; then
	if [ "$PIGZ_INSTALLED" = "YES" ]; then
		ZIPPER="pigz"
	else
		ZIPPER="gzip"
	fi
	echo $(date +%H:%M) "Compressing fasta and fastq files..."
	find "${ANALYSIS_DIR}" -type f -name '*.fasta' -exec ${ZIPPER} "{}" \;
	find "${ANALYSIS_DIR}" -type f -name '*.fastq' -exec ${ZIPPER} "{}" \;
	echo $(date +%H:%M) "Cleanup performed."
else
	echo $(date +%H:%M) "Cleanup not performed."
fi

FINISH_TIME=$(date +%Y%m%d_%H%M)
curl -sS http://textbelt.com/text -d number=$PHONE_NUMBER -d message="Pipeline finished! Started $START_TIME Finished $FINISH_TIME"

echo -e '\nAll finished!'
