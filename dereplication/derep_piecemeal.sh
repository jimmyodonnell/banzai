#!/usr/bin/env bash

# this script requires four arguments:

# 1:
# The path to a text file (e.g. nosingle.txt), where each line contains a name, DNA sequence, and semicolon-separated sequence headers that contain the library "lib" and primer "tag" index information.
# Example: >DUP_1;2;AATTCGAATT;3:21403:26840:9032_lib_J_tag_GCGCTC; 3:21404:1193:9047_lib_J_tag_GCGCTC
infile="${1}"

# 2:
# The file path to the sequencing metadata, to grab "lib" and "tag" info
SEQUENCING_METADATA="${2}"

# 3:
# The name of the column that contains the LIBRARY information for each sample
LIBRARY_COLUMN_NAME="${3}"

# 4:
# The name of the column that contains the TAG information for each sample
TAG_COLUMN_NAME="${4}"


# Get the number of available cores. This allows for parallelization without overloading.
n_cores=$(getconf _NPROCESSORS_ONLN)

# assign a path for the output (a table of counts of each duplicate sequence in each unique combination of library and primer ("tag") indexes
duplicate_table="${infile%/*}"/duplicate_table.csv

# make a directory to store the temporary duplicate files
temp_dir="${infile%/*}"/dup_temp
mkdir "${temp_dir}"


# read the library index column
LIB_COL=$(awk -F',' -v LIB_COL_NAME=$LIBRARY_COLUMN_NAME '{for (i=1;i<=NF;i++) if($i == LIB_COL_NAME) print i; exit}' $SEQUENCING_METADATA)
LIBS=$(awk -F',' -v LIBCOL=$LIB_COL 'NR>1 {print $LIBCOL}' $SEQUENCING_METADATA | sort | uniq)
echo $LIBS

# read the primer tag index column
TAG_COL=$(awk -F',' -v TAG_COL_NAME=$TAG_COLUMN_NAME '{for (i=1;i<=NF;i++) if($i == TAG_COL_NAME) print i; exit}' $SEQUENCING_METADATA)
TAGS=$(awk -F',' -v TAGCOL=$TAG_COL 'NR>1 {print $TAGCOL}' $SEQUENCING_METADATA | sort | uniq)
echo $TAGS

LIB_TAG_MOD=$(awk -F',' -v LIBCOL=$LIB_COL -v TAGCOL=$TAG_COL 'NR>1 {print "lib_" $LIBCOL "_tag_" $TAGCOL}' $SEQUENCING_METADATA | sort | uniq)

echo $LIB_TAG_MOD

# say the process has started
echo $(date +%H:%M) "Consolidating identical sequences per sample... (awk)"

# for each of the library + tag combinations,
# for libtag in $LIB_TAG_MOD; do

# set a file prefix for the batches of samples
sample_batch_prefix="${infile%/*}"/sample_batch_

# split the sample identifiers (lib + tag combination) into batches of no more than the number of available cores
echo $LIB_TAG_MOD | tr ' ' '\n' | split -l "${n_cores}" - "${sample_batch_prefix}"

################################################################################
# NEW
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
		awk 'BEGIN {print "'$sample'" ; FS ="'${sample}'" } { print NF -1 }' "${infile}" > "${temp_dir}"/"${sample}".dup
		
		) &
	
	done
	
	wait
	
done

# exit

# rm "${sample_batch_prefix}"*
################################################################################



################################################################################
# OLD
# for each of the library indexes,
# for CURRENT_LIB in $LIBS; do
	# Isolate this process to be put in the background
# 	(	
	# and for each of the primer tags in each of the libraries,
# 	for TAG_SEQ in $TAGS; do	
		# smash together the currently processed library and primer tag index as they are expected to occur in the file (e.g. "lib_A_tag_AATGAT")
# 		LIB_TAG=lib_"${CURRENT_LIB}_tag_${TAG_SEQ}"		
		# say that it's being processed
# 		echo $(date +%H:%M) "Processing" "${LIB_TAG}""..."
		# write an output file called *.dup, start by printing the lib/tag being processed, then print a count the occurrences of the current lib/tag on each line of the input file
# 		awk 'BEGIN {print "'$LIB_TAG'" ; FS ="'${LIB_TAG}'" } { print NF -1 }' "${infile}" > "${temp_dir}"/"${LIB_TAG}".dup
# 		( awk 'BEGIN {print "'$libtag'" ; FS ="'${libtag}'" } { print NF -1 }' "${infile}" > "${temp_dir}"/"${libtag}".dup ) &
# 	done
	# and put that process in the background	
# 	) &
# done
# wait
################################################################################


# write a file of names of each of the duplicates:
awk -F';' 'BEGIN {print "sample"} {print $1}' $infile | sed 's/>//' > ${temp_dir}/dupnames

# first, count the number of duplicate files:
n_files=$(find "${temp_dir}" -type f -name '*.dup*' | wc -l)
max_files=$(ulimit -n)

# find all of the *.dup files, and paste them together, and put a name at the beginning of each line.
# find "${temp_dir}" -type f -name '*.dup' -exec paste -d, {} \+ | awk '{ print "DUP_" NR-1 "," $0 }' > "${duplicate_table}"
# after increasing ulimit (maximum number of files allowed to be open at once), this took "user	1m21.943s"
# on 20150806, I got the error: paste: dup_temp/lib_J_tag_ATGCAG.dup: Too many open files

# this pastes vertically (all one column)
# for file in "${temp_dir}"/*; do
# # 	echo $file
# 	paste -d, $file >> "${duplicate_table}"
# done

# This fails
# paste -d, "${temp_dir}"/*.dup >> "${duplicate_table}"


# this will paste row by row... takes 48s on a set of 300 files (samples) each containing 630023 lines (duplicates)
paste -s -d, "${temp_dir}"/dupnames "${temp_dir}"/*.dup > "${duplicate_table}"

# this will do columns; it takes a very long time.
# for file in "${temp_dir}"/*; do cat final.dup | paste - $file >temp; cp temp final.dup; done; rm temp

# delete all of the '.dup' files
# rm "${infile%/*}"/*.dup

# say that you're finished.
echo $(date +%H:%M) "Identical sequences consolidated in file ""${duplicate_table}"

# Write fasta file in order to blast sequences
# echo $(date +%H:%M) "Writing fasta file of duplicate sequences"
awk -F';' '{ print $1 ";size=" $2 ";\n" $3 }' "${infile}" > "${infile%/*}"/duplicates.fasta

