#!/usr/bin/env bash

# usage:
# bash "/path/to/blast_script.sh" "/path/to/input/query.fasta"

# Suggestions below are based on tests run by Ryan Kelly and Jimmy O'Donnell

# QUERY
# a fasta file, read as the first argument
input_fasta="${1}"

# DATABASE
# for full nt on UW CEG server: blast_db="/local/blast-local-db/nt"
blast_db="/local/blast-local-db/nt"

# OUTPUT FORMAT
# suggested outputs: XML (5) or uncommented tabular (6) or commented tabular (7)
# "5" # XML
# "6 qseqid sallseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" # ***readable by MEGAN***
# "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" # default from http://www.ncbi.nlm.nih.gov/books/NBK279675/: 
# "7 qseqid sseqid pident staxids sscinames scomnames sblastnames"
output_format="6 qseqid sallseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle"

# IDENTITY
# percent identity suggestions: 97, 98, 99

# NUMBER OF MATCHES
# suggested: 200, 500

# CULLING_LIMIT
# suggestion changed from 5 to 20 (20150805) because the lower (and default) number can produce odd results when there are several species with similar high scores.

# WORD SIZE
# larger word sizes yield substantial speedups. RPK suggests 30.

# E VALUE
# No suggestions as of 20150819

# Automatically detect the time and set it to make a unique filename
start_time=$(date +%Y%m%d_%H%M)

# Automatically detect and set the number of cores
n_cores=$(getconf _NPROCESSORS_ONLN)

blastn \
	-db "${blast_db}" \
	-query "${input_fasta}" \
	-perc_identity 95 \
	-word_size 30 \
	-evalue 1e-20 \
	-max_target_seqs 500 \
	-culling_limit 20 \
	-outfmt  "${output_format}" \
	-out blasted_"${start_time}".tsv \
	-num_threads "${n_cores}"



################################################################################
# GRAVEYARD
################################################################################

# blastn -db /local/blast-local-db/nt -query 9_OTUs.fasta -perc_identity 97 -word_size 11 -evalue 1e-20 -max_target_seqs 500 -outfmt  5 -out blasted_"${start_time}".xml -num_threads 4 &


# old Usage:
# nohup bash 'blast_script.sh file.fasta' &
# "nohup" and "&" allow you to run it on the server, close your connection, and leave it running.

