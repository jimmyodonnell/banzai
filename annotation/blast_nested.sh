#!/usr/bin/env bash

# usage:
# bash "/path/to/blast_script.sh" "/path/to/input/query.fasta"

# Suggestions below are based on tests run by Ryan Kelly and Jimmy O'Donnell


# Automatically detect the time and set it to make a unique filename
start_time=$(date +%Y%m%d_%H%M)

echo $(date +%H:%M) "Running nested BLAST script..."

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



################################################################################
#~~~~~~~~~~~~~~~~~~~~ BLAST PARAMETERS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
################################################################################

# IDENTITY
# percent identity suggestions: 97, 98, 99
nested_identities=( 100 99 98 97 95 90 85 80 )
echo "Blast will run at these identity values:" ${nested_identities[@]}

# DATABASE
# full nt on UW CEG server: blast_db="/local/blast-local-db/nt"
# full nt on NWFSC iMac: /Users/jimmy.odonnell/NCBI/databases/nt/nt
blast_db="/Users/jimmy.odonnell/NCBI/databases/nt/nt"

# OUTPUT FORMAT
# suggested outputs: XML (5) or uncommented tabular (6) or commented tabular (7)
# "5" # XML
# "6 qseqid sallseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" # ***readable by MEGAN***
# "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" # default from http://www.ncbi.nlm.nih.gov/books/NBK279675/:
# "7 qseqid sseqid pident staxids sscinames scomnames sblastnames"
output_format="6 qseqid sallseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle"

# NUMBER OF MATCHES
# suggested: 200, 500
num_matches="500"

# CULLING_LIMIT
# "If the query range of a hit is enveloped by that of at least this many higher-scoring hits, delete the hit"
# suggestion changed from 5 to 20 (20150805) because the lower (and default) number can produce odd results when there are several species with similar high scores.
culling_limit="20"

# WORD SIZE
# larger word sizes yield substantial speedups. Smaller words yield more hits.
# default = 11; minimum = 7
# RPK suggests 30.
word_size="10"

# E VALUE
# No suggestions as of 20150819
evalue="1e-20"

# Automatically detect and set the number of cores
n_cores=$(getconf _NPROCESSORS_ONLN)

################################################################################


# assemble output directory path
out_dir="${fasta_orig%/*}"/"blast_""${start_time}"
echo "Output will be stored in:"
echo "${out_dir}"
echo

# make the output directory
mkdir "${out_dir}"

# number of identities
# N_iter="${#nested_identities[@]}"

for iter in "${nested_identities[@]}"
# for iter in "${nested_identities}"
do

	# input file
	if [ $iter == "${nested_identities[0]}" ]; then
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

	echo "${iter}" "${iter}" "${iter}" "${iter}" "${iter}" "${iter}" "${iter}"
	echo "Performing blast search at" "${iter}""% identity"


	echo "blast will query file:" "($N_seq sequence(s))"
	echo "${fasta_iter}"

	# echo "input sequence IDs (fasta headers) in file:"
	# echo "${infile_seqids}"
	# echo

	# make the output file name based on the choice of format
	outfile_base="${out_dir}"/blasted_"${start_time}"
	if [[ "${output_format}" = "5" ]] ; then
		extension="xml"
	else
		extension="txt"
	fi
	outfile="${outfile_base}"_i"${iter}"."${extension}"

	echo "blastn is running at identity" $iter "... (started at $(date +%H:%M))"

	blastn \
		-db "${blast_db}" \
		-query "${fasta_iter}" \
		-perc_identity "${iter}" \
		-word_size "${word_size}" \
		-evalue "${evalue}" \
		-max_target_seqs "${num_matches}" \
		-culling_limit "${culling_limit}" \
		-outfmt  "${output_format}" \
		-out "${outfile}" \
		-num_threads "${n_cores}"


	echo "Finished blast at "${iter}"% identity at $(date +%H:%M)"
	echo "Blast results in file:"
	echo "${outfile}"
	echo
	# touch $outfile

	# EXTRACT SEQUENCES THAT HAD NO HITS
	# make file path
	no_hits="${outfile_base}"_i"${iter}".nohits

	# grab the sequence IDs from the blast outfile, remove duplicates, compare to the seqIDs in the input fasta file, write to new file
	awk '{ print $1 }' "${outfile}" | uniq | comm -31 - "${infile_seqids}" > "${no_hits}"
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
	fasta_out="${outfile_base}"_i"${iter}"_nohits.fasta
	echo "Sequences with no blast hit can be found in fasta file:"
	echo "${fasta_out}"
	echo
	# touch $fasta_out

	grep -f "${no_hits}" "${fasta_iter}" -A 1 | sed '/^--$/d' > "${fasta_out}"

	# rm "${infile_seqids}"

done

# exit
