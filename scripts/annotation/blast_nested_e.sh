#!/usr/bin/env bash

# bash script requiring two arguments:
# 1. path to fasta file to blast
# 2. quoted and space separated e-values at which to perform nested blast.
# e.g. : "4.430273e-52 3.077510e-48 2.137807e-44 1.485037e-40 1.031588e-36 7.165977e-33 4.977879e-29 3.457907e-25 2.402052e-21 1.668597e-17 1.159098e-13"
# note: is probably NOT impervious to duplicate or unordered values

# Standard usage:
# bash "/path/to/blast_script.sh" "/path/to/input/query.fasta" "0 4.430273e-52 3.077510e-48 2.137807e-44 1.485037e-40 1.031588e-36 7.165977e-33 4.977879e-29 3.457907e-25 2.402052e-21 1.668597e-17 1.159098e-13"

# Automatically detect the time and set it to make a unique filename
start_time=$(date +%Y%m%d_%H%M)
dir_name=blast_"${start_time}"

# Optional upload to drive:
# Optional third argument: the ID (long string) of the Google Drive folder to which you'd like to upload results.
# EG: 0B_rFWkh8Szupd1J6US16aW1SMDQ
# EG: 0B_rFWkh8SzupWFNrNDVKVHZNRE0 (temp)
if [ "${3}" ] ; then
	gdrive_parent="${3}"
	# create a subdirectory for the output of this run, and grab its google drive ID
	results_dirID_gdrive=$(drive folder --title "${dir_name}" --parent "${gdrive_parent}" | awk 'NR==1 { print $2 } ')
	echo "Google Drive parent directory ID read"
	if hash "drive" 2>/dev/null; then
	# if command -v "${i}" >/dev/null 2>&1; then
		echo 'Found Google Drive command line program in' "drive" 'in' $( which "drive" )
		use_googledrive="YES"
		echo 'This confirms the blast job began on '"$(date)"'.
		Results will be uploaded to this directory as they are produced.
		The quality of the results may be enhanced by leaving out some cookies and milk for the programming elf who wrote this script, while it runs and you sleep.' | \
		drive upload --stdin --title 'blast_initiated.txt' --parent "${results_dirID_gdrive}"

	else
		echo 'Google Drive parent directory ID given, but executable not found'
		echo 'Google Drive command line executable available here:'
		echo 'https://github.com/prasmussen/gdrive'
		echo 'Must be called "drive" and be stored in a directory in PATH'
	fi

else
	echo "No Google Drive parent directory ID given"
fi





# Suggestions below are based on tests run by Ryan Kelly and Jimmy O'Donnell



# echo $(date) "Running nested BLAST script..."

# QUERY
# a fasta file, read as the first argument
fasta_orig="${1}"
# fasta_orig="/Users/jimmy.odonnell/Desktop/temp.fasta"

# check if input file exists:
if [[ -s "${fasta_orig}" ]] ; then
	echo "Original fasta input:"
	echo "${fasta_orig}"
	echo
else
	echo
	echo 'ERROR! Could not find input fasta file. You specified the file path:'
	echo
	echo "${fasta_orig}"
	echo
	echo 'That file is empty or does not exist. Aborting script.'
	exit
fi

# assemble output directory path
out_dir="${fasta_orig%/*}"/"blast_""${start_time}"
echo "Output will be stored in this directory:"
echo "${out_dir}"
echo

# make the output directory
mkdir "${out_dir}"

# make the logfile
LOGFILE="${out_dir}"/logfile.log
exec > >(tee "${LOGFILE}") 2>&1


################################################################################
echo ~~~~~~~~~~~~~~~~~~~~ BLAST PARAMETERS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
################################################################################

# IDENTITY
# percent identity suggestions: 97, 98, 99 (lower is fine with nested e value)
identity="0"
echo "Identity value:" ${identity}

# DATABASE
# full nt on UW CEG server: blast_db="/local/blast-local-db/nt"
# full nt on NWFSC iMac: /Users/jimmy.odonnell/NCBI/databases/nt/nt
blast_db="/Users/jimmy.odonnell/NCBI/databases/nt/nt"

if [[ -d "${blast_db%/*}" ]]; then
  echo "Database:" "${blast_db}"
else
	echo "${blast_db%/*}" "is not a directory on this computer. Aborting script."
	exit
fi


# NUMBER OF MATCHES
# suggested: 200, 500
num_matches="1000"
echo "Maximum matches per sequence:" "${num_matches}"

# CULLING_LIMIT
# "If the query range of a hit is enveloped by that of at least this many higher-scoring hits, delete the hit"
# suggestion changed from 5 to 20 (20150805) because the lower (and default) number can produce odd results when there are several species with similar high scores.
culling_limit="100"
echo "Max high quality hits before culling lower quality:" "${culling_limit}"

# WORD SIZE
# larger word sizes yield substantial speedups. Smaller words yield more hits.
# default = 11; minimum = 7
# RPK suggests 30.
word_size="7"
echo "Word size:" "${word_size}"

# E VALUE
# No suggestions as of 20150819
# nested e values: "1e-100"
# grab from argument 2:
if [[ -n "${2}" ]]; then
	# read the argument into an array
	arg_evalue_array=( $(echo "${2}") )
	# sort the array by decreasing numeric value and grab only uniq values
	# nested_evalues=($(printf '%s\n' "${arg_evalue_array[@]}" | sort -nr | uniq ))
	# don't sort - e values are in exponential notation, which causes problems
	nested_evalues=($(printf '%s\n' "${arg_evalue_array[@]}" ))
	echo "Nested e-values read from command line argument."
else
	nested_evalues=( 5e-52 5e-48 5e-44 5e-40 5e-36 5e-33 5e-29 5e-25 5e-21 5e-17 5e-13 )
	echo "Nested e-values argument not found; setting to default."
fi
echo "E-Values:" "${nested_evalues[@]}"

# Automatically detect and set the number of cores
n_cores=$(getconf _NPROCESSORS_ONLN)
echo "Number of threads:" "${n_cores}"

# OUTPUT FORMAT
# suggested outputs: XML (5) or uncommented tabular (6) or commented tabular (7)
# "5" # XML
# "6 qseqid sallseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" # ***readable by MEGAN***
# "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" # default from http://www.ncbi.nlm.nih.gov/books/NBK279675/:
# "7 qseqid sseqid pident staxids sscinames scomnames sblastnames"
output_format="6 qseqid sallseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids stitle"
echo "Output format:" "${output_format}"

################################################################################

# number of identities
# N_iter="${#nested_evalues[@]}"

for iter in "${nested_evalues[@]}"
# for iter in "${nested_evalues}"
do

	# input file
	if [ "${iter}" == "${nested_evalues[0]}" ]; then
		fasta_iter="${fasta_orig}"
	else
		fasta_iter="${nohits_fasta}"
	fi

	# count the number of input sequences
	N_seq=$( grep '^>' -c "${fasta_iter}" )

	infile_seqids="${fasta_iter%.*}".names

	# note that leading '>' and trailing ';' are removed
	awk '/^>/ { print }' $fasta_iter | sed -e 's/^>//g' -e 's/;$//g' > "${infile_seqids}"
	# alt: awk '/^>/ { print }' $fasta_iter | sort | uniq

	# echo "input sequence IDs (fasta headers) in file:"
	# echo "${infile_seqids}"
	# echo

	# make the output file name based on the choice of format
	hits_base="${out_dir}"/blasted_"${start_time}"
	if [[ "${output_format}" = "5" ]] ; then
		extension="xml"
	else
		extension="txt"
	fi
	hits="${hits_base}"_e"${iter}"."${extension}"

	echo "${iter}" "${iter}" "${iter}" "${iter}" "${iter}" "${iter}" "${iter}"
	echo "blast will query file:" "($N_seq sequence(s))"
	echo "${fasta_iter}"
	echo $(date) "blastn is running at" $iter" e-value..."
	time1=$(date -u +"%s")

	blastn \
		-db "${blast_db}" \
		-query "${fasta_iter}" \
		-perc_identity "${identity}" \
		-word_size "${word_size}" \
		-evalue "${iter}" \
		-max_target_seqs "${num_matches}" \
		-culling_limit "${culling_limit}" \
		-outfmt  "${output_format}" \
		-out "${hits}" \
		-num_threads "${n_cores}"


	echo $(date) "blastn finished at "${iter}" e-value."
	time2=$(date -u +"%s")
	diff=$(($time2-$time1))
	echo "blast search took $(($diff / 3600)) hours, $((($diff / 60) % 60)) minutes and $(($diff % 60)) seconds."

	echo "Blast results in file:"
	echo "${hits}"
	echo
	# touch $hits

	# EXTRACT SEQUENCES THAT HAD NO HITS
	# make file path
	no_hits="${hits_base}"_e"${iter}".nohits

	# grab the sequence IDs from the blast hits, remove duplicates, compare to the seqIDs in the input fasta file, write to new file
	grep -f <(awk '{ print $1 }' "${hits}" | uniq ) "${infile_seqids}" -v > "${no_hits}"

	# alt 0: (results in 'sequence leakage'(!) probably due to improper sorting before the comm command) awk '{ print $1 }' "${hits}" | uniq | comm -31 - "${infile_seqids}" > "${no_hits}"
	# alt 1: awk '{ print $1 }' $blast_out | sort | uniq
	# alt 2: cut -d '    ' -f 1

	if [ "${use_googledrive}" = "YES" ]; then
		drive upload --file "${hits}" --parent "${results_dirID_gdrive}"
	fi

	if [[ -s "${no_hits}" ]]; then
		N_nohits=$( wc -l "${no_hits}" | awk '{ print $1 }' )
		echo "${N_nohits}" "Sequences had no hit in database. Sequence IDs with no blast hit can be found in file:"
		echo "${no_hits}"
		echo
	else
		echo "All sequences had hits at e-value" "${iter}"
		echo "Exiting nested blast analysis."
		break
	fi



	# Write fasta file for next blast:
	nohits_fasta="${hits_base}"_e"${iter}"_nohits.fasta
	echo "Sequences with no blast hit can be found in fasta file:"
	echo "${nohits_fasta}"
	echo
	# touch $nohits_fasta

	grep -f "${no_hits}" "${fasta_iter}" -A 1 | sed '/^--$/d' > "${nohits_fasta}"
	if [ "${use_googledrive}" = "YES" ]; then
		drive upload --file "${nohits_fasta}" --parent "${results_dirID_gdrive}"
	fi

	# rm "${infile_seqids}"

done

echo $(date) "Nested e-value blast completed."

if [ "${use_googledrive}" = "YES" ]; then
	drive upload --file "${LOGFILE}" --parent "${results_dirID_gdrive}"
fi

# exit
