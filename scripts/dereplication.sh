#!/usr/bin/env bash

################################################################################
# CONSOLIDATE IDENTICAL SEQUENCES (DEREPLICATION)
################################################################################
echo $(date +%Y-%m-%d\ %H:%M) "Identifying identical sequences... (python)"
derep_output="${DEREP_INPUT}".derep
python "$SCRIPT_DIR/scripts/dereplication/dereplicate_fasta.py" "${DEREP_INPUT}"

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
echo $(date +%Y-%m-%d\ %H:%M) "Counting duplicates per identical sequence and excluding singletons if specified... (awk)"

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
echo $ID_COMBO | tr ' ' '\n' | split -l "${n_cores}" - "${sample_batch_prefix}"


# for each of the batches of files
for batch in "${sample_batch_prefix}"* ; do

	# echo processing "${batch##*/}"

	# 	current_batch=$( cat "${batch}" ) # this reads whitespace rather than newline

	for sample in $( cat "$batch" ) ; do

		# say that it's being processed
		echo $(date +%Y-%m-%d\ %H:%M) "Processing" "${sample}""..."

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
echo $(date +%Y-%m-%d\ %H:%M) "Identical sequences consolidated in file:"
echo "${duplicate_table}"
echo

# Write fasta file in order to blast sequences
echo $(date +%Y-%m-%d\ %H:%M) "Writing fasta file of duplicate sequences"
awk -F';' '{ print $1 ";size=" $2 ";\n" $3 }' "${dup_counts}" > "${duplicate_fasta}"
echo

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
