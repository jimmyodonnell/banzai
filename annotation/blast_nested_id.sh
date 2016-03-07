#!/usr/bin/env bash

# bash script requiring two arguments:
# 1. path to fasta file to blast
# 2. quoted and space separated identity values at which to perform nested blast
# e.g. : "100 99 98 97 95 90 85 80" (default if no value given.)

# note: should be impervious to duplicate or unordered values

# usage:
# bash "/path/to/blast_script.sh" "/path/to/input/query.fasta" "100 99 98 97 95 90 85 80"

# Suggestions below are based on tests run by Ryan Kelly and Jimmy O'Donnell


# Automatically detect the time and set it to make a unique filename
start_time=$(date +%Y%m%d_%H%M)

# echo $(date +%H:%M) "Running nested BLAST script..."

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
LOGFILE="${out_dir}"/logfile.txt
exec > >(tee "${LOGFILE}") 2>&1


################################################################################
echo ~~~~~~~~~~~~~~~~~~~~ BLAST PARAMETERS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
################################################################################

# IDENTITY
# percent identity suggestions: 97, 98, 99
# grab from argument 2:
if [[ -n "${2}" ]]; then
	# read the argument into an array
	arg_id_arr=( $(echo "${2}") )
	# sort the array by decreasing numeric value and grab only uniq values
	nested_identities=($(printf '%s\n' "${arg_id_arr[@]}" | sort -nr | uniq ))
	echo "Nested identities read from command line argument."
else
	nested_identities=( 100 99 98 97 95 90 85 80 )
	echo "Nested identities argument not found; setting to default."
fi
echo "Identity values:" ${nested_identities[@]}

# DATABASE
# full nt on UW CEG server: blast_db="/local/blast-local-db/nt"
# full nt on NWFSC iMac: /Users/jimmy.odonnell/NCBI/databases/nt/nt
blast_db="/Users/jimmy.odonnell/NCBI/databases/nt/nt"
echo "Database:" "${blast_db}"

# NUMBER OF MATCHES
# suggested: 200, 500
num_matches="1000"
echo "Maximum matches per sequence:" "${num_matches}"

# CULLING_LIMIT
# "If the query range of a hit is enveloped by that of at least this many higher-scoring hits, delete the hit"
# suggestion changed from 5 to 20 (20150805) because the lower (and default) number can produce odd results when there are several species with similar high scores.
culling_limit="100"
echo "Max high quality hits before culling lower quality:" "${blast_db}"

# WORD SIZE
# larger word sizes yield substantial speedups. Smaller words yield more hits.
# default = 11; minimum = 7
# RPK suggests 30.
word_size="9"
echo "Word size:" "${word_size}"

# E VALUE
# No suggestions as of 20150819
# nested e values: "1e-100"
evalue="1e-20"
echo "E-Value:" "${evalue}"

# Automatically detect and set the number of cores
n_cores=$(getconf _NPROCESSORS_ONLN)
echo "Number of threads:" "${n_cores}"

# OUTPUT FORMAT
# suggested outputs: XML (5) or uncommented tabular (6) or commented tabular (7)
# "5" # XML
# "6 qseqid sallseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" # ***readable by MEGAN***
# "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" # default from http://www.ncbi.nlm.nih.gov/books/NBK279675/:
# "7 qseqid sseqid pident staxids sscinames scomnames sblastnames"
output_format="6 qseqid sallseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle"
echo "Output format:" "${output_format}"

################################################################################

# number of identities
# N_iter="${#nested_identities[@]}"

for iter in "${nested_identities[@]}"
# for iter in "${nested_identities}"
do

	# input file
	if [ "${iter}" == "${nested_identities[0]}" ]; then
		fasta_iter="${fasta_orig}"
	else
		fasta_iter="${fasta_out}"
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
	hits="${hits_base}"_i"${iter}"."${extension}"

	echo "${iter}" "${iter}" "${iter}" "${iter}" "${iter}" "${iter}" "${iter}"
	echo "blast will query file:" "($N_seq sequence(s))"
	echo "${fasta_iter}"
	echo $(date +%H:%M) "blastn is running at" $iter"% identity..."

	blastn \
		-db "${blast_db}" \
		-query "${fasta_iter}" \
		-perc_identity "${iter}" \
		-word_size "${word_size}" \
		-evalue "${evalue}" \
		-max_target_seqs "${num_matches}" \
		-culling_limit "${culling_limit}" \
		-outfmt  "${output_format}" \
		-out "${hits}" \
		-num_threads "${n_cores}"


	echo $(date +%H:%M) "blastn finished at "${iter}"% identity."
	echo "Blast results in file:"
	echo "${hits}"
	echo
	# touch $hits

	# EXTRACT SEQUENCES THAT HAD NO HITS
	# make file path
	no_hits="${hits_base}"_i"${iter}".nohits

	# grab the sequence IDs from the blast hits, remove duplicates, compare to the seqIDs in the input fasta file, write to new file
	grep -f <(awk '{ print $1 }' "${hits}" | uniq ) "${infile_seqids}" -v > "${no_hits}"

	# alt 0: (results in 'sequence leakage'(!) probably due to improper sorting before the comm command) awk '{ print $1 }' "${hits}" | uniq | comm -31 - "${infile_seqids}" > "${no_hits}"
	# alt 1: awk '{ print $1 }' $blast_out | sort | uniq
	# alt 2: cut -d '    ' -f 1

	if [[ -s "${no_hits}" ]]; then
		N_nohits=$( wc -l "${no_hits}" | awk '{ print $1 }' )
		echo "${N_nohits}" "Sequences had no hit in database. Sequence IDs with no blast hit can be found in file:"
		echo "${no_hits}"
		echo
	else
		echo "All sequences had hits at identity" "${iter}"
		echo "Exiting nested blast analysis."
		break
	fi


	# Write fasta file for next blast:
	fasta_out="${hits_base}"_i"${iter}"_nohits.fasta
	echo "Sequences with no blast hit can be found in fasta file:"
	echo "${fasta_out}"
	echo
	# touch $fasta_out

	grep -f "${no_hits}" "${fasta_iter}" -A 1 | sed '/^--$/d' > "${fasta_out}"

	# rm "${infile_seqids}"

done

echo $(date +%H:%M) "Nested percent identity blast completed."

# exit
