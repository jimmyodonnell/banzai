#!/usr/bin/env bash

# Pipeline for analysis of MULTIPLEXED Illumina data, a la Jimmy

echo
echo
echo -e '\t' "\x20\xf0\x9f\x8f\x84" " "  "\xc2\xa1" BANZAI !
echo
echo


################################################################################
# CHECK FOR RAW DATA
################################################################################

# Define a variable called START_TIME
START_TIME=$(date +%Y%m%d_%H%M)
START_TIME_SEC=$(date +%Y%m%d_%H%M%S)
time1=$(date -u +"%s")

# Find the directory this script lives in, so it can find its friends.
SCRIPT_DIR="$(dirname "$0")"

# Specify the parameter file, check for test mode
if [[ "${1}" == "test" ]]; then
	param_file="${SCRIPT_DIR}"/banzai_params.sh
else
	param_file="${1}"
fi

# check if param file exists:
if [[ -s "${param_file}" ]] ; then
	echo "Reading analysis parameters from:"
	echo "${param_file}"
	source "${param_file}"
	echo
else
	echo
	echo 'ERROR! Could not find analysis parameter file. You specified the file path:'
	echo
	echo "${param_file}"
	echo
	echo 'That file is empty or does not exist. Aborting script.'
	exit
fi

# override singleton call to ensure test conformity
if [[ "${1}" == "test" ]]; then
	remove_singletons="YES"
fi

# check if sequencing metadata exists
if [[ -s "${SEQUENCING_METADATA}" ]] ; then
	echo "Reading sequencing metadata from:"
	echo "${SEQUENCING_METADATA}"
	echo
else
	echo
	echo 'ERROR! Could not find sequencing metadata file. You specified the file path:'
	echo
	echo "${SEQUENCING_METADATA}"
	echo
	echo 'That file is empty or does not exist. Aborting script.'
	exit
fi

# check for correct newline characters (CRLF will break things)
source "${SCRIPT_DIR}"/scripts/newline_fix.sh "${SEQUENCING_METADATA}"
if [[ -s "${NEWLINES_FIXED}" ]]; then
	SEQUENCING_METADATA="${NEWLINES_FIXED}"
fi

################################################################################
# CHECK FOR DEPENDENCIES
################################################################################
dependencies=($( echo pear cutadapt vsearch swarm seqtk python blastn R ))
source "${SCRIPT_DIR}"/scripts/dependency_check.sh "${dependencies[@]}"

# Load functions
for file in "${SCRIPT_DIR}"/scripts/functions/* ; do
	source "${file}"
done

# Specify compression utility
if hash pigz 2>/dev/null; then
	ZIPPER="pigz"
	echo "pigz installation found"
	echo
else
	ZIPPER="gzip"
	echo "pigz installation not found; using gzip"
	echo
fi

# Detect number of cores on machine; set variable
n_cores=$(getconf _NPROCESSORS_ONLN)
if [ $n_cores -gt 1 ]; then
	echo "$n_cores cores detected."
	echo
else
	n_cores=1
	echo "Multiple cores not detected."
	echo
fi

# make an analysis directory with starting time timestamp
OUTPUT_DIR="${OUTPUT_DIRECTORY}"/banzai_out_"${START_TIME}"
if [[ -d "${OUTPUT_DIR}" ]]; then
	OUTPUT_DIR="${OUTPUT_DIRECTORY}"/banzai_out_"${START_TIME_SEC}"
	if [[ -d "${OUTPUT_DIR}" ]]; then
		echo "Output directory already exists!"
		echo "${OUTPUT_DIR}"
		echo "Aborting script."
		exit
	fi
fi
mkdir "${OUTPUT_DIR}"

# Write a log file of output from this script (everything that prints to terminal)
LOGFILE="${OUTPUT_DIR}"/logfile.txt
exec > >(tee "${LOGFILE}") 2>&1

echo $(date +%Y-%m-%d\ %H:%M) "Analysis started at ""${START_TIME}"
echo "Output is located in:"
echo "${OUTPUT_DIR}"
echo

# Copy these files into that directory as a verifiable log you can refer back to.
cp "${SCRIPT_DIR}"/banzai.sh "${OUTPUT_DIR}"/analysis_script.txt
cp "${param_file}" "${OUTPUT_DIR}"/analysis_parameters.txt



################################################################################
# READ METADATA
################################################################################
# report metadata dimensions
METADATA_DIM=($( awk -F, 'END{print NR, NF}' "${SEQUENCING_METADATA}" ))
echo "Metadata has" "${METADATA_DIM[0]}" "rows and" "${METADATA_DIM[1]}" "columns including header."
N_SAMPLES=$( echo "${METADATA_DIM[0]}" - 1 | bc )
echo "Expecting" "${N_SAMPLES}" "samples total."
echo

# Filnames
COLNUM_FILE1=$( get_colnum "${COLNAME_FILE1}" "${SEQUENCING_METADATA}")
COLNUM_FILE2=$( get_colnum "${COLNAME_FILE2}" "${SEQUENCING_METADATA}")

# Library names
COLNUM_ID1=$( get_colnum "${COLNAME_ID1_NAME}" "${SEQUENCING_METADATA}")
COLNUM_ID1_SEQ=$( get_colnum "${COLNAME_ID1_SEQ}" "${SEQUENCING_METADATA}")

if [[ "${SECONDARY_INDEX}" == "YES" ]]; then
	# Secondary multiplex index sequences
	COLNUM_ID2=$( get_colnum "${COLNAME_ID2_SEQ}" "${SEQUENCING_METADATA}")

	# Secondary index sequence positions
	COLNUM_ID2_START=$( get_colnum "${COLNAME_ID2_START}" "${SEQUENCING_METADATA}")

fi

# Sample names
COLNUM_SAMPLE=$( get_colnum "${COLNAME_SAMPLE_ID}" "${SEQUENCING_METADATA}")

# Primers
COLNUM_PRIMER1=$( get_colnum "${COLNAME_PRIMER1}" "${SEQUENCING_METADATA}")
COLNUM_PRIMER2=$( get_colnum "${COLNAME_PRIMER2}" "${SEQUENCING_METADATA}")

################################################################################
# CHECK FILES
################################################################################
FILE1=($(awk -F',' -v COLNUM=$COLNUM_FILE1 \
  'NR>1 {  print $COLNUM }' $SEQUENCING_METADATA |\
  sort | uniq))

FILE2=($(awk -F',' -v COLNUM=$COLNUM_FILE2 \
  'NR>1 {print $COLNUM}' $SEQUENCING_METADATA |\
  sort | uniq ))

NFILE1="${#FILE1[@]}"
NFILE2="${#FILE2[@]}"
if [ "${NFILE1}" != "${NFILE2}" ]; then
	echo "ERROR: Whoa! different number of forward and reverse files"
fi

if [[ -n "${FILE1}" && -n "${FILE2}" ]]; then
  echo 'Files read from metadata columns' "${COLNUM_FILE1}" 'and' "${COLNUM_FILE2}"
  echo 'File names:'
	for (( i=0; i < "${NFILE1}"; ++i)); do
		printf '%s\t%s\n' "${FILE1[i]}" "${FILE2[i]}"
	done
	echo
else
  echo 'ERROR:' 'At least one file is not valid'
  echo 'Looked in metadata columns' "${COLNUM_FILE1}" 'and' "${COLNUM_FILE2}"
  echo 'Aborting script'
  exit
fi

################################################################################
# LOAD SECONDARY INDEXES
################################################################################
if [[ "${SECONDARY_INDEX}" == "YES" ]]; then

	ID2S=($(awk -F',' -v COLNUM=$COLNUM_ID2 \
	  'NR>1 {  print $COLNUM }' $SEQUENCING_METADATA |\
	  sort | uniq))
	N_index_sequences="${#ID2S}"
	ID2_LENGTH=${#ID2S[0]}
	ID2_START=($(awk -F',' -v COLNUM=$COLNUM_ID2_START \
	  'NR>1 {  print $COLNUM }' $SEQUENCING_METADATA |\
	  sort | uniq))

	# check if number of indexes is greater than one:
	if [[ "${N_index_sequences}" -gt 1 ]]; then
		echo "Secondary indexes read from sequencing metadata (""${N_index_sequences}"" total)"
		echo
	else
	  echo
	  echo 'ERROR:' "${N_index_sequences}" 'index sequences found. There should probably be more than 1.'
	  echo
	  echo 'Aborting script.'
		exit
	fi

fi

################################################################################
# Read in primers
################################################################################
PRIMER1=($(awk -F',' -v COLNUM=$COLNUM_PRIMER1 \
  'NR > 1 { print $COLNUM }' $SEQUENCING_METADATA |\
  sort | uniq ))

PRIMER2=($(awk -F',' -v COLNUM=$COLNUM_PRIMER2 \
  'NR > 1 { print $COLNUM }' $SEQUENCING_METADATA |\
  sort | uniq ))

if [[ -n "${PRIMER1}" && -n "${PRIMER2}" ]]; then
  echo 'Primers read from metadata columns' "${COLNUM_PRIMER1}" 'and' "${COLNUM_PRIMER2}"
  echo 'Primer sequences:' "${PRIMER1}" "${PRIMER2}"
	echo
else
  echo 'ERROR:' 'At least one primer is not valid'
  echo 'Looked in metadata columns' "${COLNUM_PRIMER1}" 'and' "${COLNUM_PRIMER2}"
  echo 'Aborting script'
  exit
fi


################################################################################
# Calculate the expected size of the region of interest
################################################################################
EXTRA_SEQ=${ID2S[0]}${ID2S[0]}$PRIMER1$PRIMER2
LENGTH_ROI=$(( $LENGTH_FRAG - ${#EXTRA_SEQ} ))
LENGTH_ROI_HALF=$(( $LENGTH_ROI / 2 ))


# Unique samples are given by combining the primary and secondary indexes
ID_COMBO=$( awk -F',' -v COLNUM1=$COLNUM_ID1 -v COLNUM2=$COLNUM_ID2 \
'NR>1 {
  print ";ID1=" $COLNUM1 ";ID2=" $COLNUM2
}' "${SEQUENCING_METADATA}" )

SAMPLE_NAMES=($(awk -F',' -v COLNUM=$COLNUM_SAMPLE \
  'NR>1 { print $COLNUM }' "${SEQUENCING_METADATA}" ))

ID1_ALL=($(awk -F',' -v COLNUM=$COLNUM_ID1 \
  'NR>1 { print $COLNUM }' "${SEQUENCING_METADATA}" ))

if [[ "${SECONDARY_INDEX}" == "YES" ]]; then
  ID2_ALL=($(awk -F',' -v COLNUM=$COLNUM_ID2 \
    'NR>1 { print $COLNUM }' "${SEQUENCING_METADATA}" ))
fi

ID2_ALL_RC=($( for i in "${ID2_ALL[@]}"; do revcom $i; done))

# write file for translating demultiplexed output to samples
SAMPLE_TRANS_FILE="${OUTPUT_DIR}"/sample_trans.tmp
for (( i=0; i < "${#ID2_ALL[@]}"; i++ )); do
  printf "ID1=%s;ID2A=%s;ID2B=%s\tID1=%s;ID2=%s\tsample=%s\n" \
	"${ID1_ALL[i]}" "${ID2_ALL[i]}" "${ID2_ALL_RC[i]}" \
	"${ID1_ALL[i]}" "${ID2_ALL[i]}" \
	"${SAMPLE_NAMES[i]}" >> "${SAMPLE_TRANS_FILE}"
done

# create a file to store index efficiency data
INDEX_COUNT="${OUTPUT_DIR}"/index_count.txt
echo "index1 index2 left_side right_side" >> "${INDEX_COUNT}"

# create a directory to store concatenated output
CONCAT_DIR="${OUTPUT_DIR}"/all_lib
mkdir "${CONCAT_DIR}"
CONCAT_FILE="${CONCAT_DIR}"/demult_concat.fasta

################################################################################
# BEGIN LOOP TO PERFORM LIBRARY-LEVEL ACTIONS
################################################################################

for (( i=0; i < "${#FILE1[@]}"; i++ )); do

	# Identify the forward and reverse fastq files.
	CURRENT_FILE[0]="${FILE1[i]}"
	CURRENT_FILE[1]="${FILE2[i]}"

	# PEAR v0.9.6 does not correctly merge .gz files.
	# Look through files and decompress if necessary.
	for j in 0 1 ; do
	  FILEPATH=$( find "${PARENT_DIR}" -name "${CURRENT_FILE[j]}"*  )
		if file "${FILEPATH}" | grep -q gzip ; then
			echo $(date +%Y-%m-%d\ %H:%M) "decompressing "${FILEPATH}""
			"${ZIPPER}" -d "${FILEPATH}"
			READ[j]=$( find "${PARENT_DIR}" -name "${CURRENT_FILE[j]}"*  )
		else
			READ[j]="${FILEPATH}"
		fi
  done

	READ1="${READ[0]}"
	READ2="${READ[1]}"

	for j in "${READ1}" "${READ2}"; do
		if [[ ! -s "${j}" ]] ; then
		    echo 'Could not find the raw reads file! Looked here:'
				echo "${j}"
		    echo 'Aborting script.'
		    exit
		fi
  done

  CURRENT_ID1_SEQ=$( awk -F, '
	/'"${CURRENT_FILE[0]}"'/ { print $'"${COLNUM_ID1_SEQ}"'; }' "${SEQUENCING_METADATA}" |\
	sort | uniq )

  CURRENT_ID1_NAME=$( awk -F, '
	/'"${CURRENT_FILE[0]}"'/ { print $'"${COLNUM_ID1}"'; }' "${SEQUENCING_METADATA}" |\
	sort | uniq )

	ID1_OUTPUT_DIR="${OUTPUT_DIR}"/${CURRENT_ID1_NAME}
	mkdir "${ID1_OUTPUT_DIR}"

	##############################################################################
	# MERGE PAIRED-END READS AND QUALITY FILTER (PEAR)
	##############################################################################

	LENGTH_READ=$( head -n 100000 "${READ1}" | awk '{print length($0);}' |\
	  sort -nr | uniq | head -n 1 )

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
		echo
	else
		echo $(date +%Y-%m-%d\ %H:%M) "Merging reads in library" "${CURRENT_ID1_NAME}""..."
		MERGED_READS_PREFIX="${ID1_OUTPUT_DIR}"/1_merged
		MERGED_READS="${ID1_OUTPUT_DIR}"/1_merged.assembled.fastq
		if [[ "${USE_PEAR_DEFAULTS}" == "NO" ]]; then
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
		else
			pear \
			--forward-fastq "${READ1}" \
			--reverse-fastq "${READ2}" \
			--output "${MERGED_READS_PREFIX}"
    fi

		# check pear output:
		if [[ ! -s "${MERGED_READS}" ]] ; then
		    echo 'ERROR: No reads were merged.'
		    echo 'Aborting analysis of this library, but will move on to next one.'
				continue
		fi

		if [ "${HOARD}" = "NO" ]; then
			rm "${ID1_OUTPUT_DIR}"/1_merged.discarded.fastq
			rm "${ID1_OUTPUT_DIR}"/1_merged.unassembled.forward.fastq
			rm "${ID1_OUTPUT_DIR}"/1_merged.unassembled.reverse.fastq
		fi

		echo

	fi

	################################################################################
	# EXPECTED ERROR FILTERING (vsearch)
	################################################################################
	# FILTER READS (This is the last step that uses quality scores, so convert to fasta)
	if [ "${Perform_Expected_Error_Filter}" = "YES" ]; then
		echo $(date +%Y-%m-%d\ %H:%M) "Filtering merged reads..."
		FILTERED_OUTPUT="${ID1_OUTPUT_DIR}"/2_filtered.fasta
		vsearch \
			--fastq_filter "${MERGED_READS}" \
			--fastq_maxee "${Max_Expected_Errors}" \
			--fastaout "${FILTERED_OUTPUT}" \
			--fasta_width 0

    echo
	else
		# Convert merged reads fastq to fasta
		echo  $(date +%Y-%m-%d\ %H:%M) "converting fastq to fasta..."
		FILTERED_OUTPUT="${MERGED_READS%.*}".fasta
		seqtk seq -A "${MERGED_READS}" > "${FILTERED_OUTPUT}"
		echo
	fi

	if [ "${HOARD}" = "NO" ]; then
		rm "${MERGED_READS}"
	fi

	# Compress merged reads
  echo $(date +%Y-%m-%d\ %H:%M) "Compressing PEAR output..."
  find "${ID1_OUTPUT_DIR}" -type f -name '*.fastq' -exec ${ZIPPER} "{}" \;
  echo $(date +%Y-%m-%d\ %H:%M) "PEAR output compressed."
	echo


	FILTERED_RENAMED="${FILTERED_OUTPUT%.*}"_renamed.fasta
	if [ "${RENAME_READS}" = "YES" ]; then
		echo $(date +%Y-%m-%d\ %H:%M) "Renaming reads in library" "${CURRENT_ID1_NAME}""..."

		awk -F'[: ]' '{
				if ( /^>/ )
					print ">"$4":"$5":"$6":"$7";ID1='${CURRENT_ID1_NAME}'";
				else
					print $0
		}' "${FILTERED_OUTPUT}" > "${FILTERED_RENAMED}"

		mv "${FILTERED_RENAMED}" "${FILTERED_OUTPUT}"

		echo $(date +%Y-%m-%d\ %H:%M) "Reads renamed"
		echo

	else

		awk '{
				if ( /^>/ )
					print $0";ID1='${CURRENT_ID1_NAME}'";
				else
					print $0
		}' "${FILTERED_OUTPUT}" > "${FILTERED_RENAMED}"

		mv "${FILTERED_RENAMED}" "${FILTERED_OUTPUT}"

		echo "Reads not renamed"
		echo

	fi


	################################################################################
	# HOMOPOLYMERS (grep, awk)
	################################################################################
	if [ "${REMOVE_HOMOPOLYMERS}" = "YES" ]; then
		echo $(date +%Y-%m-%d\ %H:%M) "Removing homopolymers..."
		NOHOMO_FASTA="${ID1_OUTPUT_DIR}"/3_no_homopolymers.fasta
		source "${SCRIPT_DIR}"/scripts/homopolymer_rm.sh "${FILTERED_OUTPUT}" "${HOMOPOLYMER_MAX}" "${NOHOMO_FASTA}"
		DEMULTIPLEX_INPUT="${NOHOMO_FASTA}"
	else
		echo "Homopolymers not removed."
		DEMULTIPLEX_INPUT="${FILTERED_OUTPUT}"
		echo
	fi

	################################################################################
	# DEMULTIPLEXING (awk)
	################################################################################
	if [[ "${SECONDARY_INDEX}" == "YES" ]]; then
	  CURRENT_ID1_DEMULT="${ID1_OUTPUT_DIR}"/demult.fasta
		source "${SCRIPT_DIR}"/beta/demulti.sh -i "${DEMULTIPLEX_INPUT}" \
		  -s "${ID2_START}" -l "${ID2_LENGTH}" >> "${CURRENT_ID1_DEMULT}"
		if [ "${HOARD}" = "NO" ]; then
			rm "${DEMULTIPLEX_INPUT}"
		fi

		echo $(date +%Y-%m-%d\ %H:%M) "Concatenating fasta files..."
		cat "${ID1_OUTPUT_DIR}"/demult.fasta >> "${CONCAT_FILE}"
		echo
	else
		cat "${DEMULTIPLEX_INPUT}" >> "${CONCAT_FILE}"
  fi

	################################################################################
	# ID1 SPECIFIC PRIMER REMOVAL
	################################################################################
	# ID1_PRIMER_REM="${ID1_OUTPUT_DIR}"/demult.fasta
	# cat "${ID1_OUTPUT_DIR}"/demultiplexed/*/2_notags.fasta >> "${ID1_PRIMER_REM}"
	ID1_PRIMER_REMOVAL_OUT="${ID1_OUTPUT_DIR}"/no_primers.fasta
	# source "${SCRIPT_DIR}"/scripts/primer_removal.sh \
	#   "${CURRENT_ID1_DEMULT}" "${ID1_PRIMER_REMOVAL_OUT}" \
	#   "${PRIMER1}"  "${PRIMER2}" \
	# 	"${PRIMER_MISMATCH_PROPORTION}"  "${LENGTH_ROI_HALF}"

	if [ "${HOARD}" = "NO" ]; then
		rm -rf "${ID1_OUTPUT_DIR}"
	fi

done

################################################################################
# END LOOP TO PERFORM LIBRARY-LEVEL ACTIONS
################################################################################
PRIMER_REM_ID1_ALL="${CONCAT_DIR}"/no_primer_byid1.fasta
# cat "${OUTPUT_DIR}"/*/demult.fasta >> "${PRIMER_REM_ID1_ALL}"

################################################################################
# PRIMER REMOVAL
################################################################################
PRIMER_REMOVAL_OUT="${CONCAT_DIR}"/no_primers.fasta
source "${SCRIPT_DIR}"/scripts/primer_removal.sh \
  "${CONCAT_FILE}" "${PRIMER_REMOVAL_OUT}" \
  "${PRIMER1}"  "${PRIMER2}" \
	"${PRIMER_MISMATCH_PROPORTION}"  "${LENGTH_ROI_HALF}"

# check for cutadapt/primer removal success.
if [[ ! -s "${PRIMER_REMOVAL_OUT}" ]]; then
  echo 'ERROR: cutadapt did not process reads correctly. This file is empty or absent:'
	echo "${PRIMER_REMOVAL_OUT}"
  echo 'Aborting script'
	echo
  exit
fi

################################################################################
# CONSOLIDATE IDENTICAL SEQUENCES (DEREPLICATION)
################################################################################
echo $(date +%Y-%m-%d\ %H:%M) "Identifying identical sequences... (python)"
DEREP_FASTA="${CONCAT_DIR}"/derep.fasta
DEREP_MAP="${CONCAT_DIR}"/derep.map
python "$SCRIPT_DIR/scripts/dereplication/derep_fasta.py" \
  "${PRIMER_REMOVAL_OUT}" 'ID1=' "${DEREP_FASTA}" "${DEREP_MAP}"

# check if duplicate fasta and duplicate table exist. (Might need to check size)
if [[ ! -s "${DEREP_FASTA}" ]] ; then
    echo 'There was a problem generating the duplicate fasta file. It is empty or absent.'
    echo 'The remainder of the script depends on this file.'
    echo 'Aborting script.'
    exit
fi
# check if duplicate fasta and duplicate table exist. (Might need to check size)
if [[ ! -s "${DEREP_MAP}" ]] ; then
    echo 'There was a problem generating the dereplication map file. It is empty or absent.'
    echo 'The remainder of the script depends on this file.'
    echo 'Aborting script.'
    exit
fi

DUPLICATE_TABLE="${CONCAT_DIR}"/duplicate_table.csv

# Rscript "${SCRIPT_DIR}"/scripts/dereplication/long_to_wide.R \
#   "${DEREP_MAP}" "${DUPLICATE_TABLE}" "${remove_singletons}"

# check if duplicate table exists. (Might need to check size)
if [[ ! -s "${DEREP_MAP}" ]] ; then
    echo 'There was a problem generating the duplicate table. It is empty or absent.'
    echo 'Aborting script.'
    exit
fi

##############################################################################
# CHECK FOR CHIMERAS
##############################################################################
if [[ "${remove_chimeras}" = "YES" ]] ; then
echo $(date +%Y-%m-%d\ %H:%M) 'Looking for chimeras in duplicate fasta file using vsearch'
source "${SCRIPT_DIR}"/scripts/chimera_check.sh "${DEREP_FASTA}"
clustering_input="${chimera_free_fasta}"
echo
else
clustering_input="${DEREP_FASTA}"
fi





################################################################################
# CLUSTER OTUS
################################################################################
# Note that identical (duplicate) sequences were consolidated earlier;
# This step outputs a file (*.uc) that lists, for every sequence, which sequence it clusters with
if [ "$CLUSTER_OTUS" = "NO" ]; then
	BLAST_INPUT="${clustering_input}"
else
	case "${cluster_method}" in

	    "swarm" )

	        echo $(date +%Y-%m-%d\ %H:%M) 'Clustering sequences into OTUs using swarm'
	        source "${SCRIPT_DIR}"/scripts/OTU_clustering/cluster_swarm.sh "${clustering_input}"
					echo

	    ;;

	    "vsearch" )

	        # echo $(date +%Y-%m-%d\ %H:%M) 'Clustering sequences into OTUs using vsearch'
	        # source "${SCRIPT_DIR}"/scripts/OTU_clustering/cluster_vsearch.sh "${duplicate_fasta}"
					echo "Sorry, OTU clustering with vsearch has not been implemented yet."
					echo $(date +%Y-%m-%d\ %H:%M) 'Clustering sequences into OTUs using swarm'
	        source "${SCRIPT_DIR}"/scripts/OTU_clustering/cluster_swarm.sh "${clustering_input}"
					echo

	    ;;

	    "usearch" )

	        echo $(date +%Y-%m-%d\ %H:%M) 'Clustering sequences into OTUs using usearch'
	        source "${SCRIPT_DIR}"/scripts/OTU_clustering/cluster_usearch.sh "${clustering_input}"
					echo

	    ;;

	    * )

	        echo "${cluster_method}" 'is an invalid clustering method.'
	        echo 'Must be one of swarm, vsearch, usearch, or none.'
	        echo $(date +%Y-%m-%d\ %H:%M) 'Clustering sequences into OTUs using swarm'
	        source "${SCRIPT_DIR}"/scripts/OTU_clustering/cluster_swarm.sh "${clustering_input}"
					echo


	    ;;

	esac

	# check that dup to otu map is greater than 12 bytes
	minsize=12
	size_dup_otu_map=$(wc -c <"${dup_otu_map}")
	if [ $size_dup_otu_map -lt $minsize ]; then
	    echo 'There was an error generating the dup-to-otu map.'
			echo
	fi


	# Assign the path for the OTU to sample map file
	OTU_map="${dup_otu_map%/*}"/OTU.map

	# Convert duplicate table to OTU table using R script (arguments: 
  # (1) unique sequence ("duplicate") to sample map file, 
  # (2) unique sequence ("duplicate") to otu map file, 
  # (3) otu to sample map file path
  Rscript "$SCRIPT_DIR/scripts/dup_to_OTU.R" "${DEREP_MAP}" "${dup_otu_map}" "${OTU_map}"

	# check if OTU table and OTU fasta exist (and/or are of size gt 1?)
	if [[ ! -s "${OTU_fasta}" ]] ; then
	    echo 'There was a problem generating the OTU fasta file. It is empty or absent.'
	    echo 'Aborting script.'
	    exit
	fi
	if [[ ! -s "${OTU_map}" ]] ; then
	    echo 'There was a problem generating the OTU table. It is empty or absent.'
	    echo 'Aborting script.'
	    exit
	fi

fi


################################################################################
# CLEAN UP
################################################################################
if [ "$PERFORM_CLEANUP" = "YES" ]; then
	echo $(date +%Y-%m-%d\ %H:%M) "Compressing fasta, fastq, and xml files..."
	find "${OUTPUT_DIR}" -type f -name '*.fasta' -exec ${ZIPPER} "{}" \;
	find "${OUTPUT_DIR}" -type f -name '*.fastq' -exec ${ZIPPER} "{}" \;
	find "${OUTPUT_DIR}" -type f -name '*.xml' -exec ${ZIPPER} "{}" \;
	echo $(date +%Y-%m-%d\ %H:%M) "Cleanup performed."
else
	echo $(date +%Y-%m-%d\ %H:%M) "Cleanup not performed."
fi

FINISH_TIME=$(date +%Y%m%d_%H%M)
time2=$(date -u +"%s")
timediff=$(($time2-$time1))
echo "banzai took $(($timediff / 3600)) hours, $((($timediff / 60) % 60)) minutes and $(($timediff % 60)) seconds."

mail -e
if [[ "$?" < 1 ]]; then
	echo 'Pipeline finished! Started at' $START_TIME 'and finished at' $FINISH_TIME | mail -s "banzai is finished" "${EMAIL_ADDRESS}"
fi

################################################################################
# WRITE SUMMARY
################################################################################
# TODO add this; probably there is a problem caused by reading a logfile while writing to it.
# SUMMARY_FILE="${OUTPUT_DIR}"/summary.txt
# echo "Writing summary file..."
# source "${SCRIPT_DIR}"/scripts/summarize.sh "${LOGFILE}" > "${SUMMARY_FILE}"
# echo "Summary written to:"
# echo "${SUMMARY_FILE}"
# echo

################################################################################
# EXIT
################################################################################
echo -e '\n'$(date +%Y-%m-%d\ %H:%M)'\tAll finished! Why not treat yourself to a...\n'
echo
echo -e '\t~~~ MAI TAI ~~~'
echo -e '\t2 oz\taged rum'
echo -e '\t0.75 oz\tfresh squeezed lime juice'
echo -e '\t0.5 oz\torgeat'
echo -e '\t0.5 oz\ttriple sec'
echo -e '\t0.25 oz\tsimple syrup'
echo -e '\tShake, strain, and enjoy!' '\xf0\x9f\x8d\xb9\x0a''\n'

if [[ "${1}" == "test" ]]; then
  source 	"${SCRIPT_DIR}"/test/test_uniqmap.sh
  test_uniqmap "${PARENT_DIR}"/derep_expected.map "${OUTPUT_DIR}"
	if [[ "$?" < 1 ]]; then
		echo "banzai passes test"
	fi
fi
