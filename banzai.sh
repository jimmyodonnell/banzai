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

# Find the directory this script lives in, so it can find its friends.
SCRIPT_DIR="$(dirname "$0")"

# Read in the parameter file (was source "$SCRIPT_DIR/banzai_params.sh"; now argument 1)
param_file="${1}"
source "${param_file}"

# check if param file exists:
if [[ -s "${param_file}" ]] ; then
	echo "Reading analysis parameters from:"
	echo "${param_file}"
else
	echo
	echo 'ERROR! Could not find analysis parameter file. You specified the file path:'
	echo
	echo "${param_file}"
	echo
	echo 'That file is empty or does not exist. Aborting script.'
	exit
fi


# check if sequencing metadata exists
if [[ -s "${SEQUENCING_METADATA}" ]] ; then
	echo "Reading sequencing metadata from:"
	echo "${SEQUENCING_METADATA}"
else
	echo
	echo 'ERROR! Could not find sequencing metadata file. You specified the file path:'
	echo
	echo "${SEQUENCING_METADATA}"
	echo
	echo 'That file is empty or does not exist. Aborting script.'
	exit
fi

################################################################################
# CHECK FOR DEPENDENCIES
################################################################################
dependencies=($( echo pear cutadapt vsearch swarm seqtk python blastn R ))
echo 'Checking for dependencies:' "${dependencies[@]}"
for i in "${dependencies[@]}"; do
	if hash "${i}" 2>/dev/null; then
	# if command -v "${i}" >/dev/null 2>&1; then
		echo 'Found program' "${i}" 'in' $( which "${i}" )
	else
		echo 'ERROR: A program on which this script depends was not found:' "${i}"
		echo 'Aborting script.'
		exit
	fi
done

# Specify compression utility
if hash pigz 2>/dev/null; then
	ZIPPER="pigz"
	echo "pigz installation found"
else
	ZIPPER="gzip"
	echo "pigz installation not found; using gzip"
fi

# Detect number of cores on machine; set variable
n_cores=$(getconf _NPROCESSORS_ONLN)
if [ $n_cores -gt 1 ]; then
	echo "$n_cores cores detected."
else
	n_cores=1
	echo "Multiple cores not detected."
fi

# make an analysis directory with starting time timestamp
ANALYSIS_DIR="${ANALYSIS_DIRECTORY}"/Analysis_"${START_TIME}"
mkdir "${ANALYSIS_DIR}"

# Write a log file of output from this script (everything that prints to terminal)
LOGFILE="${ANALYSIS_DIR}"/logfile.txt
exec > >(tee "${LOGFILE}") 2>&1

echo "Analysis started at ""${START_TIME}" " and is located in ""${ANALYSIS_DIR}"

# Copy these files into that directory as a verifiable log you can refer back to.
cp "${SCRIPT_DIR}"/banzai.sh "${ANALYSIS_DIR}"/analysis_script.txt
cp "${param_file}" "${ANALYSIS_DIR}"/analysis_parameters.txt



################################################################################
# LOAD MULTIPLEX TAGS
################################################################################
TAG_COL=$(awk -F',' -v TAG_COL_NAME=$TAG_COLUMN_NAME '{for (i=1;i<=NF;i++) if($i == TAG_COL_NAME) print i; exit}' $SEQUENCING_METADATA)
TAGS=$(awk -F',' -v TAGCOL=$TAG_COL 'NR>1 {print $TAGCOL}' $SEQUENCING_METADATA | sort | uniq)
N_index_sequences=$(echo $TAGS | awk '{print NF}')

# check if number of tags is greater than one:
if [[ "${N_index_sequences}" -gt 1 ]]; then
	echo "Multiplex tags read from sequencing metadata (""${N_index_sequences}"") total"
else
  echo
  echo 'ERROR:' "${N_index_sequences}" 'index sequences found. There should probably be more than 1.'
  echo
  echo 'Aborting script.'
	exit
fi

declare -a TAGS_ARRAY=($TAGS)


################################################################################
# Read in primers and create reverse complements.
################################################################################
PRIMER1_COLNUM=$(awk -F',' -v PRIMER1_COL=$PRIMER_1_COLUMN_NAME '{for (i=1;i<=NF;i++) if($i == PRIMER1_COL) print i; exit}' $SEQUENCING_METADATA)
PRIMER2_COLNUM=$(awk -F',' -v PRIMER2_COL=$PRIMER_2_COLUMN_NAME '{for (i=1;i<=NF;i++) if($i == PRIMER2_COL) print i; exit}' $SEQUENCING_METADATA)
PRIMER1=$(awk -F',' -v PRIMER1_COL=$PRIMER1_COLNUM 'NR==2 {print $PRIMER1_COL}' $SEQUENCING_METADATA)
PRIMER2=$(awk -F',' -v PRIMER2_COL=$PRIMER2_COLNUM 'NR==2 {print $PRIMER2_COL}' $SEQUENCING_METADATA)
echo "Primers read from sequencing metadata:" "${PRIMER1}" "${PRIMER2}"

if [[ -n "${PRIMER1}" && -n "${PRIMER2}" ]]; then
  echo 'Primers read from metadata columns' "${PRIMER1_COLNUM}" 'and' "${PRIMER2_COLNUM}"
  echo 'Primer sequences:' "${PRIMER1}" "${PRIMER2}"
else
  echo 'ERROR:' 'At least one primer is not valid'
  echo 'Looked in metadata columns' "${PRIMER1_COLNUM}" 'and' "${PRIMER2_COLNUM}"
  echo 'Aborting script'
  exit
fi

# make primer array
read -a primers_arr <<< $( echo $PRIMER1 $PRIMER2 )

# Reverse complement primers
PRIMER1RC=$( echo ${PRIMER1} | tr "[ABCDGHMNRSTUVWXYabcdghmnrstuvwxy]" "[TVGHCDKNYSAABWXRtvghcdknysaabwxr]" | rev )
PRIMER2RC=$( echo ${PRIMER2} | tr "[ABCDGHMNRSTUVWXYabcdghmnrstuvwxy]" "[TVGHCDKNYSAABWXRtvghcdknysaabwxr]" | rev )

# make primer array
read -a primersRC_arr <<< $( echo $PRIMER1RC $PRIMER2RC )


################################################################################
# Calculate the expected size of the region of interest, given the total size of fragments, and the length of primers and tags
################################################################################
EXTRA_SEQ=${TAGS_ARRAY[0]}${TAGS_ARRAY[0]}$PRIMER1$PRIMER2
LENGTH_ROI=$(( $LENGTH_FRAG - ${#EXTRA_SEQ} ))
LENGTH_ROI_HALF=$(( $LENGTH_ROI / 2 ))


################################################################################
# Find raw sequence files
################################################################################
# Look for any file with '.fastq' in the name in the parent directory
# note that this will include ANY file with fastq -- including QC reports!
# TODO make LIBRARY_DIRECTORIES an array by wrapping it in ()
LIBRARY_DIRECTORIES=$( find "$PARENT_DIR" -name '*.fastq*' -print0 | xargs -0 -n1 dirname | sort --unique )

# PEAR v0.9.6 does not correctly merge .gz files.
# Look through files and decompress if necessary.
raw_files=($( find "${PARENT_DIR}" -name '*.fastq*' ))
for myfile in "${raw_files[@]}"; do
	if [[ "${myfile}" =~ \.gz$ ]]; then
		echo $(date +%H:%M) "decompressing "${myfile}""
		"${ZIPPER}" -d "${myfile}"
	fi
done

# Count library directories and print the number found
# TODO if LIBRARY_DIRECTORIES is an array, its length is "${#LIBRARY_DIRECTORIES[@]}"
N_library_dir=$(echo $LIBRARY_DIRECTORIES | awk '{print NF}')
echo "${N_library_dir}"" library directories found:"

# Show the libraries that were found:
# TODO for i in "${LIBRARY_DIRECTORIES[@]}"; do echo "${i##*/}" ; done
for i in $LIBRARY_DIRECTORIES; do echo "${i##*/}" ; done

# Assign it to a variable for comparison
LIBS_FROM_DIRECTORIES=$(for i in $LIBRARY_DIRECTORIES; do echo "${i##*/}" ; done)

# Read library names from file or sequencing metadata
if [ "${READ_LIB_FROM_SEQUENCING_METADATA}" = "YES" ]; then
	LIB_COL=$(awk -F',' -v LIB_COL_NAME=$LIBRARY_COLUMN_NAME '{for (i=1;i<=NF;i++) if($i == LIB_COL_NAME) print i; exit}' $SEQUENCING_METADATA)
	LIBS=$(awk -F',' -v LIBCOL=$LIB_COL 'NR>1 {print $LIBCOL}' $SEQUENCING_METADATA | sort | uniq)
	N_libs=$(echo $LIBS | awk '{print NF}')
	echo "Library names read from sequencing metadata (""${N_libs}"") total"
	echo "${LIBS}"
else
	LIBS=$(tr '\n' ' ' < "${LIB_FILE}" )
	N_libs=$(echo $LIBS | awk '{print NF}')
	echo "Library names read from lib file (""${LIBS}"") total"
fi
# make library names into an array
# TODO LIBS_ARRAY is never used
# declare -a LIBS_ARRAY=($LIBS)

# Check that library names are the same in the metadata and file system
if [ "$LIBS_FROM_DIRECTORIES" != "$LIBS" ]; then
	echo "Warning: Library directories and library names in metadata are NOT the same. Something will probably go wrong later..."
else
	echo "Library directories and library names in metadata are the same - great job."
fi


# Unique samples are given by combining the library and tags
# TODO originally contained sort | uniq; this is unnecessary I think
LIB_TAG_MOD=$( awk -F',' -v LIBCOL=$LIB_COL -v TAGCOL=$TAG_COL 'NR>1 { print "lib_" $LIBCOL "_tag_" $TAGCOL }' $SEQUENCING_METADATA | sort | uniq )

# create a file to store tag efficiency data
TAG_COUNT="${ANALYSIS_DIR}"/tag_count.txt
echo "library tag left_tagged right_tagged" >> "${TAG_COUNT}"

################################################################################
# BEGIN LOOP TO PERFORM LIBRARY-LEVEL ACTIONS
################################################################################

for CURRENT_LIB in $LIBRARY_DIRECTORIES; do

	# Identify the forward and reverse fastq files.
	READS=($(find "${CURRENT_LIB}" -name '*.fastq*'))
	READ1="${READS[0]}"
	READ2="${READS[1]}"
	# READ1=$(find "${CURRENT_LIB}" -name '*R1*fastq*')
	# READ2=$(find "${CURRENT_LIB}" -name '*R2*fastq*')

	LIB_OUTPUT_DIR="${ANALYSIS_DIR}"/${CURRENT_LIB##*/}
	mkdir "${LIB_OUTPUT_DIR}"

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
		MERGED_READS="${PEAR_OUTPUT}/$(basename $LIB_OUTPUT_DIR)"/1_merged.assembled.fastq.gz #RPK edited 20160720
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

	################################################################################
	# EXPECTED ERROR FILTERING (vsearch)
	################################################################################
	# FILTER READS (This is the last step that uses quality scores, so convert to fasta)
	if [ "${Perform_Expected_Error_Filter}" = "YES" ]; then
		echo $(date +%H:%M) "Filtering merged reads..."
		FILTERED_OUTPUT="${LIB_OUTPUT_DIR}"/2_filtered.fasta
		vsearch \
			--fastq_filter "${MERGED_READS}" \
			--fastq_maxee "${Max_Expected_Errors}" \
			--fastaout "${FILTERED_OUTPUT}" \
			--fasta_width 0
	else
		# Convert merged reads fastq to fasta
		echo  $(date +%H:%M) "converting fastq to fasta..."
				FILTERED_OUTPUT="${LIB_OUTPUT_DIR%.*}"/1_merged.assembled.fasta  #FILTERED_OUTPUT="${MERGED_READS%.*}".fasta  #edited by RPK 20160720
		seqtk seq -A "${MERGED_READS}" > "${FILTERED_OUTPUT}"
	fi

	# Compress merged reads
  echo $(date +%H:%M) "Compressing PEAR output..."
  find "${LIB_OUTPUT_DIR}" -type f -name '*.fastq' -exec ${ZIPPER} "{}" \;
  echo $(date +%H:%M) "PEAR output compressed."





	if [ "${RENAME_READS}" = "YES" ]; then
		echo $(date +%H:%M) "Renaming reads in library" "${CURRENT_LIB##*/}""..."
		# TODO remove whitespace from sequence labels?
		# sed 's/ /_/'

		# original (usearch7):	sed -E "s/ (1|2):N:0:[0-9]/_"${CURRENT_LIB##*/}"_/" "${FILTERED_OUTPUT}" > "${CURRENT_LIB}"/tmp.fasta
		# update for usearch8, which without warning removes any part of the sequence ID following a space.
		# holy shit this ads the _"${CURRENT_LIB##*/}"_ to EVERY line
		# sed -E "s/$/_"${CURRENT_LIB##*/}"_/" "${FILTERED_OUTPUT}" > "${CURRENT_LIB}"/tmp.fasta
		# sed -E "s/>([a-zA-Z0-9-]*:){4}/>/" "${CURRENT_LIB}"/tmp.fasta > "${FILTERED_OUTPUT%.*}"_renamed.fasta
		# rm "${CURRENT_LIB}"/tmp.fasta

		# updated 20150521; one step solution using awk; removes anything after the first space!
		FILTERED_RENAMED="${FILTERED_OUTPUT%.*}"_renamed.fasta
		awk -F'[: ]' '{
				if ( /^>/ )
					print ">"$4":"$5":"$6":"$7"_lib_'${CURRENT_LIB##*/}'_";
				else
					print $0
				}' "${FILTERED_OUTPUT}" > "${FILTERED_RENAMED}"

		FILTERED_OUTPUT="${FILTERED_RENAMED}"

		# rm "${FILTERED_OUTPUT}"

	else
		echo "Reads not renamed"
	fi


	################################################################################
	# HOMOPOLYMERS (grep, awk)
	################################################################################
	if [ "${REMOVE_HOMOPOLYMERS}" = "YES" ]; then
		echo $(date +%H:%M) "Removing homopolymers..."
		HomoLineNo="${CURRENT_LIB}"/homopolymer_line_numbers.txt
		grep -E -i -B 1 -n "(A|T|C|G)\1{$HOMOPOLYMER_MAX,}" "${FILTERED_OUTPUT}" | \
			cut -f1 -d: | \
			cut -f1 -d- | \
			sed '/^$/d' > "${HomoLineNo}"
		if [ -s "${HomoLineNo}" ]; then
			DEMULTIPLEX_INPUT="${CURRENT_LIB}"/3_no_homopolymers.fasta
			awk 'NR==FNR{l[$0];next;} !(FNR in l)' "${HomoLineNo}" "${FILTERED_OUTPUT}" > "${DEMULTIPLEX_INPUT}"
			awk 'NR==FNR{l[$0];next;} (FNR in l)' "${HomoLineNo}" "${FILTERED_OUTPUT}" > "${CURRENT_LIB}"/homopolymeric_reads.fasta
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
	# 20150522 changed {0,9} to {3} to eliminate flexibility (that could result in a read being assigned to >1 sample)
	# awk '/^.{0,9}'"$TAG_SEQ"'/{if (a && a !~ /^.{0,9}'"$TAG_SEQ"'/) print a; print} {a=$0}' "${DEMULTIPLEX_INPUT}" > "${TAG_DIR}"/1_tagL_present.fasta ) &

	echo $(date +%H:%M) "Demultiplexing: removing tags and adding to sequence ID in library" "${CURRENT_LIB##*/}""..."
	for TAG_SEQ in $TAGS; do
	(	TAG_DIR="${DEMULTIPLEXED_DIR}"/tag_"${TAG_SEQ}"
		mkdir "${TAG_DIR}"
		demult_file_L="${TAG_DIR}"/1_tagL_removed.fasta
	  demult_file_R="${TAG_DIR}"/2_notags.fasta

		# Left side tag
		awk 'gsub(/^.{3}'"$TAG_SEQ"'/,"") {
			if (a && a !~ /^.{3}'"$TAG_SEQ"'/)
				print a;
			print
		} {a=$0}' "${DEMULTIPLEX_INPUT}" > "${demult_file_L}"

		# Right side tag
		TAG_RC=$( echo ${TAG_SEQ} | tr "[ATGCatgcNn]" "[TACGtacgNn]" | rev )
		awk 'gsub(/'"$TAG_RC"'.{3}$/,"") {
			if (a && a !~ /'"$TAG_RC"'.{3}$/)
				print a "tag_""'"$TAG_SEQ"'";
			print
		} {a = $0}' "${demult_file_L}" > "${demult_file_R}"

		echo "${CURRENT_LIB##*/}" "${TAG_SEQ}" $(wc -l "${demult_file_L}" | \
			awk '{ print ($1/2) }') $(wc -l "${demult_file_R}" | \
			awk '{ print ($1/2)}') >> "${TAG_COUNT}" ) &

	done

	wait

	# echo $(date +%H:%M) "Demultiplexing: removing right tag and adding tag sequence to sequence ID in library" "${CURRENT_LIB##*/}""..."
	# for TAG_SEQ in $TAGS; do
	# (	TAG_DIR="${LIB_OUTPUT_DIR}"/demultiplexed/tag_"${TAG_SEQ}"
	# 	demult_file_R="${TAG_DIR}"/2_notags.fasta
	# 	TAG_RC=$( echo ${TAG_SEQ} | tr "[ATGCatgcNn]" "[TACGtacgNn]" | rev )
	# 	# 20150522 changed {0,9} to {3} to eliminate flexibility (that could result in a read being assigned to >1 sample)
	# 	awk 'gsub(/'"$TAG_RC"'.{3}$/,"") {if (a && a !~ /'"$TAG_RC"'.{3}$/) print a "tag_""'"$TAG_SEQ"'"; print } {a = $0}' "${TAG_DIR}"/1_tagL_removed.fasta > "${demult_file_R}" ) &
	# done
	#
	# wait

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

	# TODO MOVE THE VARIABLE ASSIGNMENT TO TOP; MOVE MKDIR TO TOP OF CONCAT IF LOOP
	echo $(date +%H:%M) "Concatenating fasta files..."
	CONCAT_DIR="$ANALYSIS_DIR"/all_lib
	mkdir "${CONCAT_DIR}"
	CONCAT_FILE="${CONCAT_DIR}"/1_demult_concat.fasta

	# TODO could move this into above loop after demultiplexing?
	for CURRENT_LIB in $LIBRARY_DIRECTORIES; do

		LIB_OUTPUT_DIR="${ANALYSIS_DIR}"/${CURRENT_LIB##*/}

		for TAG_SEQ in $TAGS; do
			cat "${LIB_OUTPUT_DIR}"/demultiplexed/tag_"${TAG_SEQ}"/2_notags.fasta >> "${CONCAT_FILE}"
		done

		echo $(date +%H:%M) "Compressing fasta files..."
		find "${LIB_OUTPUT_DIR}" -type f -name '*.fasta' -exec ${ZIPPER} "{}" \;
		echo $(date +%H:%M) "fasta files compressed."

	done

	################################################################################
	# Count the occurrences of '_tag_' + the 6 characters following it in the concatenated file
	################################################################################
  # TODO !!! This will fail if there are underscores in the library names !!!
	# an attempt at making this robust to underscores
	# grep -E -o '_lib_.+?(?=_tag)_tag_.{6}' "${CONCAT_DIR}"/1_demult_concat.fasta | sed 's/_lib_//;s/_tag_/ /' | sort | uniq -c | sort -nr > "${CONCAT_DIR}"/1_demult_concat.fasta.tags

	echo $(date +%H:%M) "Counting reads associated with each sample index (primer tag)..."
	grep -E -o '_lib_[^_]*_tag_.{6}' "${CONCAT_DIR}"/1_demult_concat.fasta | sed 's/_lib_//;s/_tag_/ /' | sort | uniq -c | sort -nr > "${CONCAT_DIR}"/1_demult_concat.fasta.tags

	echo $(date +%H:%M) "Summary of sequences belonging to each sample index found in ""${CONCAT_DIR}""/1_demult_concat.fasta.tags"
	################################################################################




	################################################################################
	# PRIMER REMOVAL
	################################################################################
	# (moot for concatenated file): echo $(date +%H:%M) "Removing primers in library" "${CURRENT_LIB##*/}""..."
	# Remove PRIMER1 and PRIMER2 from the BEGINNING of the reads. NOTE cutadapt1.7+ will accept ambiguities in primers.

	# count lines in primer removal input
	echo $(date +%H:%M) "Counting sequences in primer removal input..."
	seq_N_demult_concat=$( grep -e '^>' --count "${CONCAT_FILE}" )
	echo $(date +%H:%M) "${seq_N_demult_concat}" "sequences found in file" "${CONCAT_FILE}"

	# TODO wrap in '( ) &' to force into background and allow parallel processing
	# i.e.
	# for primer in "${primers_arr[@]}"; do
	# 	( cutadapt -g ^"${primer}" -e "${PRIMER_MISMATCH_PROPORTION}" -m "${LENGTH_ROI_HALF}" --discard-untrimmed "${CONCAT_FILE}" > "${CONCAT_DIR}"/5_L"${primer}"_removed.fasta ) &
	# done
	# wait

	# remove primer 1 from left side of sequences
	primerL1_removed="${CONCAT_DIR}"/5_primerL1_removed.fasta
	( cutadapt \
		-g ^"${PRIMER1}" \
		-e "${PRIMER_MISMATCH_PROPORTION}" \
		-m "${LENGTH_ROI_HALF}" \
		--discard-untrimmed \
		"${CONCAT_FILE}" > "${primerL1_removed}" ) &

	# remove primer 2 from left side of sequences
	primerL2_removed="${CONCAT_DIR}"/5_primerL2_removed.fasta
	( cutadapt \
		-g ^"${PRIMER2}" \
		-e "${PRIMER_MISMATCH_PROPORTION}" \
		-m "${LENGTH_ROI_HALF}" \
		--discard-untrimmed \
		"${CONCAT_FILE}" > "${primerL2_removed}" ) &

	wait

	# compress left primer removal input
	echo $(date +%H:%M) "Compressing left primer removal input..."
	"${ZIPPER}" "${CONCAT_DIR}"/1_demult_concat.fasta
	echo $(date +%H:%M) "Left primer removal input compressed."

	# check for cutadapt/primer removal success.
	if [[ ! -s "${primerL1_removed}" ]]; then
	  echo 'ERROR: cutadapt did not process reads correctly. This file is empty or absent:'
		echo "${primerL1_removed}"
	  echo 'Aborting script'
	  exit
	fi
	# check for cutadapt/primer removal success.
	if [[ ! -s "${primerL2_removed}" ]]; then
	  echo 'ERROR: cutadapt did not process reads correctly. This file is empty or absent:'
		echo "${primerL2_removed}"
	  echo 'Aborting script'
	  exit
	fi

	# Remove the reverse complement of primer 1 from the right side of sequences
	primerR1_removed="${CONCAT_DIR}"/6_primerR1_removed.fasta
	( cutadapt \
		-a "${PRIMER2RC}"$ \
		-e "${PRIMER_MISMATCH_PROPORTION}" \
		-m "${LENGTH_ROI_HALF}" \
		--discard-untrimmed \
		"${primerL1_removed}" > "${primerR1_removed}" ) &

	# Remove the reverse complement of primer 2 from the right side of sequences
	primerR2_removed="${CONCAT_DIR}"/6_primerR2_removed.fasta
	( cutadapt \
		-a "${PRIMER1RC}"$ \
		-e "${PRIMER_MISMATCH_PROPORTION}" \
		-m "${LENGTH_ROI_HALF}" \
		--discard-untrimmed \
		"${primerL2_removed}" > "${primerR2_removed}" ) &

	wait

	# check for cutadapt/primer removal success.
	if [[ ! -s "${primerR1_removed}" ]]; then
		echo 'ERROR: cutadapt did not process reads correctly. This file is empty or absent:'
		echo "${primerR1_removed}"
		echo 'Aborting script'
		exit
	fi
	# check for cutadapt/primer removal success.
	if [[ ! -s "${primerR2_removed}" ]]; then
		echo 'ERROR: cutadapt did not process reads correctly. This file is empty or absent:'
		echo "${primerR2_removed}"
		echo 'Aborting script'
		exit
	fi

	# Reverse-complement the sequences in which the RC of primer 1 was found on the right side
	seqtk seq -r "${CONCAT_DIR}"/6_primerR1_removed.fasta > "${CONCAT_DIR}"/6_primerR1_removedRC.fasta

	# paste together the contents of the files that primers were removed from.
	DEREP_INPUT="${CONCAT_DIR}"/7_no_primers.fasta

	cat "${CONCAT_DIR}"/6_primerR1_removedRC.fasta "${CONCAT_DIR}"/6_primerR2_removed.fasta > "${DEREP_INPUT}"

	# check that it worked (derep input / no primers)
	if [[ ! -s "${DEREP_INPUT}" ]] ; then
	    echo 'ERROR: Input file for dereplication is empty or absent.'
	    echo 'This will cause problems for all remaining steps, so script will exit.'
	    exit
	fi



	################################################################################
	# CONSOLIDATE IDENTICAL SEQUENCES (DEREPLICATION)
	################################################################################
	echo $(date +%H:%M) "Identifying identical sequences... (python)"
	derep_output="${DEREP_INPUT}".derep
	python "$SCRIPT_DIR/dereplication/dereplicate_fasta.py" "${DEREP_INPUT}"

	# check for derep output
	if [[ ! -s "${derep_output}" ]] ; then
	    echo 'ERROR: python dereplication output is empty or absent.'
	    echo 'This will cause problems for all remaining steps, so script will exit.'
	    exit
	fi


	##############################################################################
	# COUNT SEQUENCES, REMOVE SINGLETONS
	##############################################################################
	# count the number of sequences per duplicate (print NF-1), sort them by the number of sequences per duplicate (sort -nr), and precede with a name ("DUP_X", where X is the line number), excluding singleton sequences (if NF > 2) if specified
	echo $(date +%H:%M) "Counting duplicates per identical sequence and excluding singletons if specified... (awk)"

	dup_counts="${DEREP_INPUT%/*}"/dup_counts.txt

	if [ "$remove_singletons" = "YES" ]; then
		awk -F';' '{
			if (NF > 2)
				print NF-1 ";" $0
			}' "${derep_output}" | \
		sort -nr | \
		awk -F';' '{
			print ">DUP_" NR ";" $0
		}' > "${dup_counts}"
	else
		awk -F';' '{
				print NF-1 ";" $0
			}' "${derep_output}" | \
		sort -nr | \
		awk -F';' '{
			print ">DUP_" NR ";" $0
		}' > "${dup_counts}"
	fi

	# check output
	if [[ ! -s "${dup_counts}" ]] ; then
	    echo 'There was a problem generating the dup_counts file. It is empty or absent.'
	    echo 'This will cause problems counting sequences for dereplication.'
	fi



	# COUNT OCCURRENCES PER SAMPLE (LIBRARY + TAG) PER DUPLICATE

	# assign a path for the output (a table of counts of each duplicate sequence in each unique combination of library and primer ("tag") indexes, and a fasta of all the duplicate sequences.
	duplicate_table="${DEREP_INPUT%/*}"/duplicate_table.csv
	duplicate_fasta="${DEREP_INPUT%/*}"/duplicates.fasta

	# make a directory to store the temporary duplicate files
	temp_dir="${DEREP_INPUT%/*}"/dup_temp
	mkdir "${temp_dir}"

	# set a file prefix for the batches of samples
	sample_batch_prefix="${temp_dir}"/sample_batch_

	# split the sample identifiers (lib + tag combination) into batches of no more than the number of available cores
	echo $LIB_TAG_MOD | tr ' ' '\n' | split -l "${n_cores}" - "${sample_batch_prefix}"


	# for each of the batches of files
	for batch in "${sample_batch_prefix}"* ; do

		echo processing "${batch##*/}"

		# 	current_batch=$( cat "${batch}" ) # this reads whitespace rather than newline

		for sample in $( cat "$batch" ) ; do

			# say that it's being processed
			echo $(date +%H:%M) "Processing" "${sample}""..."

			# Isolate this process to be put in the background
			(

			# write an output file called *.dup, start by printing the lib/tag being processed, then print a count the occurrences of the current lib/tag on each line of the input file
			awk 'BEGIN {print "'$sample'" ; FS ="'${sample}'" } { print NF -1 }' "${dup_counts}" > "${temp_dir}"/"${sample}".dup

			) &

		done

		wait

	done

	# write a file of names of each of the duplicates:
	dupnames="${temp_dir}"/dupnames
	awk -F';' 'BEGIN {print "sample"} {print $1}' "$dup_counts" | sed 's/>//' > "${dupnames}"

	# first, count the number of duplicate files:
	n_files=$(find "${temp_dir}" -type f -name '*.dup*' | wc -l)

	# I think this was only relevant for a different approach
	max_files=$(ulimit -n)

	# this will paste row by row... takes 48s on a set of 300 files (samples) each containing 630023 lines (duplicates)
	paste -s -d, "${dupnames}" "${temp_dir}"/*.dup > "${duplicate_table}"

	# this will do columns; it takes a very long time.
	# for file in "${temp_dir}"/*; do cat final.dup | paste - $file >temp; cp temp final.dup; done; rm temp

	# cleanup
	# rm "${infile%/*}"/*.dup
	# rm "${sample_batch_prefix}"*

	# say that you're finished.
	echo $(date +%H:%M) "Identical sequences consolidated in file ""${duplicate_table}"

	# OLD/UNSTABLE/BAD -- but generated CSV in different orientation
	# This will start as many processes as you have libraries... be careful!
	# echo $(date +%H:%M) "Consolidating identical sequences per sample... (awk)"
	# for CURRENT_LIB in $LIBRARY_DIRECTORIES; do
	# 	( for TAG_SEQ in $TAGS; do
	# 		LIB_TAG="lib_${CURRENT_LIB##*/}_tag_${TAG_SEQ}"
	# 		echo $(date +%H:%M) "Processing" "${LIB_TAG}""..."
	# 		# the output of the awk function gsub is the number of replacements, so you could use this instead... however, it appears slower?
	# 		# ( awk 'BEGIN {print "'$LIB_TAG'" } { print gsub(/"'$LIB_TAG'"/,"") }' ${DEREP_INPUT%/*}/nosingle.txt > ${DEREP_INPUT%/*}/"${LIB_TAG}".dup ) &
	# 		awk 'BEGIN {print "'$LIB_TAG'" ; FS ="'${LIB_TAG}'" } { print NF -1 }' ${DEREP_INPUT%/*}/nosingle.txt > ${DEREP_INPUT%/*}/"${LIB_TAG}".dup
	# 	done ) &
	# done
	# wait
	# # Write a csv file of the number of occurrences of each duplicate sequence per tag. (rows = sequences, cols = samples)
	# duplicate_table="${DEREP_INPUT%/*}"/dups.csv
	# find "${DEREP_INPUT%/*}" -type f -name '*.dup' -exec paste -d, {} \+ | awk '{ print "DUP_" NR-1 "," $0 }' > "${duplicate_table}"
	# # delete all of the '.dup' files
	# rm ${DEREP_INPUT%/*}/*.dup

	# Write fasta file in order to blast sequences
	echo $(date +%H:%M) "Writing fasta file of duplicate sequences"
	awk -F';' '{ print $1 ";size=" $2 ";\n" $3 }' "${dup_counts}" > "${duplicate_fasta}"

	# check if duplicate fasta and duplicate table exist. (Might need to check size)
	if [[ ! -s "${duplicate_fasta}" ]] ; then
	    echo 'There was a problem generating the duplicate fasta file. It is empty or absent.'
	    echo 'The remainder of the script, including OTU clustering, depends on this file.'
	    echo 'Aborting script.'
	    exit
	fi
	if [[ ! -s "${duplicate_table}" ]] ; then
	    echo 'There was a problem generating the duplicate table. It is empty or absent.'
	    echo 'Aborting script.'
	    exit
	fi

################################################################################





##############################################################################
# CHECK FOR CHIMERAS
##############################################################################
if [[ "${remove_chimeras}" = "YES" ]] ; then
	echo $(date +%H:%M) 'Looking for chimeras in duplicate fasta file using vsearch'
	source "${SCRIPT_DIR}"/chimera_check.sh "${duplicate_fasta}"
	clustering_input="${chimera_free_fasta}"
else
	clustering_input="${duplicate_fasta}"
fi





	################################################################################
	# CLUSTER OTUS
	################################################################################
	# Note that identical (duplicate) sequences were consolidated earlier;
	# This step outputs a file (*.uc) that lists, for every sequence, which sequence it clusters with
	if [ "$CLUSTER_OTUS" = "NO" ]; then
		BLAST_INPUT="${clustering_input}"
	else
		echo $(date +%H:%M) "Clustering OTUs..."

		case "${cluster_method}" in

		    "swarm" )

		        echo $(date +%H:%M) 'Clustering sequences into OTUs using swarm'
		        source "${SCRIPT_DIR}"/OTU_clustering/cluster_swarm.sh "${clustering_input}"

		    ;;

		    "vsearch" )

		        # echo $(date +%H:%M) 'Clustering sequences into OTUs using vsearch'
		        # source "${SCRIPT_DIR}"/OTU_clustering/cluster_vsearch.sh "${duplicate_fasta}"
						echo "Sorry, OTU clustering with vsearch has not been implemented yet."
						echo $(date +%H:%M) 'Clustering sequences into OTUs using swarm'
		        source "${SCRIPT_DIR}"/OTU_clustering/cluster_swarm.sh "${clustering_input}"

		    ;;

		    "usearch" )

		        echo $(date +%H:%M) 'Clustering sequences into OTUs using usearch'
		        source "${SCRIPT_DIR}"/OTU_clustering/cluster_usearch.sh "${clustering_input}"

		    ;;

		    * )

		        echo "${cluster_method}" 'is an invalid clustering method.'
		        echo 'Must be one of swarm, vsearch, usearch, or none.'
		        echo $(date +%H:%M) 'Clustering sequences into OTUs using swarm'
		        source "${SCRIPT_DIR}"/OTU_clustering/cluster_swarm.sh "${clustering_input}"

		    ;;

		esac

		# check that dup to otu map is greater than 12 bytes
		minsize=12
		size_dup_otu_map=$(wc -c <"${dup_otu_map}")
		if [ $size_dup_otu_map -lt $minsize ]; then
		    echo 'There was an error generating the dup-to-otu map.'
		fi


		# Assign the path for the OTU table
		# OTU_table="${dir_out}"/OTU_table.csv

		# Convert duplicate table to OTU table using R script (arguments: (1) duplicate table, (2) dup to otu table, (3) otu table path
		Rscript "$SCRIPT_DIR/dup_to_OTU_table.R" "${duplicate_table}" "${dup_otu_map}" "${OTU_table}"

		# check if OTU table and OTU fasta exist (and/or are of size gt 1?)
		if [[ ! -s "${OTU_fasta}" ]] ; then
		    echo 'There was a problem generating the OTU fasta file. It is empty or absent.'
		    echo 'Aborting script.'
		    exit
		fi
		if [[ ! -s "${OTU_table}" ]] ; then
		    echo 'There was a problem generating the OTU table. It is empty or absent.'
		    echo 'Aborting script.'
		    exit
		fi

	fi

	################################################################################
	# BLAST CLUSTERS
	################################################################################
	echo $(date +%H:%M) "BLASTing..."
	blast_output="${DEREP_INPUT%/*}"/10_BLASTed.xml
	blastn \
		-query "${BLAST_INPUT}" \
		-db "$BLAST_DB" \
		-num_threads "$n_cores" \
		-perc_identity "${PERCENT_IDENTITY}" \
		-word_size "${WORD_SIZE}" \
		-evalue "${EVALUE}" \
		-max_target_seqs "${MAXIMUM_MATCHES}" \
		-culling_limit "${CULLING}" \
		-outfmt 5 \
		-out "${blast_output}"

	# check for blast output
	if [[ ! -s "${blast_output}"  ]]; then
		echo
		echo 'BLAST failed: the output file is empty or absent.'
	    echo 'File should be:' "${blast_output}"
		echo
	fi















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
		python "$SCRIPT_DIR/dereplication/dereplicate_fasta.py" "${DEREP_INPUT}"

		# REMOVE SINGLETONS
		# usearch -sortbysize "${TAG_DIR}"/7_derep.fasta -minsize 2 -sizein -sizeout -output "${TAG_DIR}"/8_nosingle.fasta
		# COUNT DUPLICATES PER READ, REMOVE SINGLETONS
		awk -F';' '{ if (NF > 2) print NF-1 ";" $0 }' "${DEREP_INPUT}".derep | sort -nr | awk -F';' '{ print ">DUP_" NR ";" $0}' > ${DEREP_INPUT%/*}/dup_counts.txt

		# count the duplicates
		awk 'BEGIN { FS ="_tag_'${TAG_SEQ}'" } { print NF -1 }' "${DEREP_INPUT%/*}"/dup_counts.txt > ${DEREP_INPUT%/*}/"${TAG_SEQ}".dup

		# Write fasta file in order to blast sequences
		awk -F';' '{ print $1 ";size=" $2 ";\n" $3 }' ${DEREP_INPUT%/*}/dup_counts.txt > ${DEREP_INPUT%/*}/no_duplicates.fasta

		# CLUSTER SEQUENCES
		if [ "$CLUSTER_OTUS" = "NO" ]; then
			BLAST_INPUT=${DEREP_INPUT%/*}/no_duplicates.fasta
		else
			CLUSTER_RADIUS="$(( 100 - ${CLUSTERING_PERCENT} ))"
			UPARSE_OUT="${DEREP_INPUT%/*}"/OTU_uparse.txt
			usearch -cluster_otus "${DEREP_INPUT%/*}"/dup_counts.txt -otu_radius_pct "${CLUSTER_RADIUS}" -sizein -sizeout -otus "${TAG_DIR}"/9_OTUs.fasta -uparseout "${UPARSE_OUT}"
			BLAST_INPUT="${TAG_DIR}"/9_OTUs.fasta
		fi

		# BLAST CLUSTERS
		blastn -query "${BLAST_INPUT}" -db "$BLAST_DB" -num_threads "$n_cores" -perc_identity "${PERCENT_IDENTITY}" -word_size "${WORD_SIZE}" -evalue "${EVALUE}" -max_target_seqs "${MAXIMUM_MATCHES}" -outfmt 5 -out "${TAG_DIR}"/10_BLASTed.xml
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

		# check for blast output
		if [[ -s "${blast_output}"  ]]; then

		echo $(date +%H:%M) 'BLAST output found; proceeding to MEGAN.'
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
		echo "export what=CSV format=readname_taxonname separator=comma file=${DIR}/meganout_${COLLAPSE_RANK1}.csv;" >> "${MEGAN_COMMAND_FILE}"
		if [ "$PERFORM_SECONDARY_MEGAN" = "YES" ]; then
			echo "collapse rank='$COLLAPSE_RANK2';" >> "${MEGAN_COMMAND_FILE}"
			echo "update;" >> "${MEGAN_COMMAND_FILE}"
			echo "select nodes=all;" >> "${MEGAN_COMMAND_FILE}"
			echo "export what=CSV format=readname_taxonname separator=comma file=${DIR}/meganout_${COLLAPSE_RANK2}.csv;" >> "${MEGAN_COMMAND_FILE}"
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

	else
		echo
		echo 'BLAST failed: the output file is empty or absent.'
		echo 'File should be:' "${blast_output}"
		echo
	fi

done


################################################################################
# PRELIMINARY ANALYSES
################################################################################
# Once you have a final CSV file of the number of occurences of each OTU in each sample, run some preliminary analyses in R
# TODO rename preliminary to OTU analyses; move analysis script to OTU analysis directory
OUTPUT_PDF="${ANALYSIS_DIR}"/analysis_results_"${START_TIME}".pdf

echo $(date +%H:%M) "passing args to R for preliminary analysis..."
Rscript "$SCRIPT_DIR/analyses_prelim.R" "${OUTPUT_PDF}" "${OTU_table}" "${SEQUENCING_METADATA}" "${LIBRARY_COLUMN_NAME}" "${TAG_COLUMN_NAME}" "${ColumnName_SampleName}" "${ColumnName_SampleType}"


# EMPTY PDFs are 3829 bytes
minimumsize=4000
size_PDF=$(wc -c <"${OUTPUT_PDF}")
if [ "${size_PDF}" -lt "${minimumsize}" ]; then
    echo 'There was a problem generating the PDF.'
else
	REMOTE_PDF="${OUTPUT_PDF_DIR}"/analysis_results_"${START_TIME}".pdf
	cp "${OUTPUT_PDF}" "${REMOTE_PDF}"
fi


if [ "$PERFORM_CLEANUP" = "YES" ]; then
	echo $(date +%H:%M) "Compressing fasta, fastq, and xml files..."
	find "${ANALYSIS_DIR}" -type f -name '*.fasta' -exec ${ZIPPER} "{}" \;
	find "${ANALYSIS_DIR}" -type f -name '*.fastq' -exec ${ZIPPER} "{}" \;
	find "${ANALYSIS_DIR}" -type f -name '*.xml' -exec ${ZIPPER} "{}" \;
	echo $(date +%H:%M) "Cleanup performed."
else
	echo $(date +%H:%M) "Cleanup not performed."
fi

FINISH_TIME=$(date +%Y%m%d_%H%M)

echo 'Pipeline finished! Started at' $START_TIME 'and finished at' $FINISH_TIME | mail -s "banzai is finished" "${EMAIL_ADDRESS}"

echo -e '\n'$(date +%H:%M)'\tAll finished! Why not treat yourself to a...\n'
echo
echo -e '\t~~~ MAI TAI ~~~'
echo -e '\t2 oz\taged rum'
echo -e '\t0.75 oz\tfresh squeezed lime juice'
echo -e '\t0.5 oz\torgeat'
echo -e '\t0.5 oz\ttriple sec'
echo -e '\t0.25 oz\tsimple syrup'
echo -e '\tShake, strain, and enjoy!' '\xf0\x9f\x8d\xb9\x0a''\n'
