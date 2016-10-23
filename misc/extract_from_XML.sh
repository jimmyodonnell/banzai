#!/usr/bin/env bash

# In what directory should I look?
parent_dir="/Users/threeprime/Desktop/Reanalysis_20141031_1601"

xml_files=($(find "${parent_dir}" -name "*.xml" 2>/dev/null))
xml_files=($(find "${parent_dir}" -name "*.xml"))

xml_files=$( find "${parent_dir}" -type f -name '*.xml' )

declare -a xml_array=($xml_files)

# what is the string of characters I should look for?
string_to_match='No hits found'

# how many lines above it is the sequence identifier? (seems like its always 15)
lines_above=15

# what do ALL of the sequence IDs start with?
# if sequence names haven't been changed, and they're all from the same run, this should be fine because it's the ID of the machine the sequence came from
seq_ID_prefix='M00'
# otherwise, you could use the XML sequence name identifier: <Iteration_query-def>

# What should I call the files containing the names and the sequences of the sequences with no hits?
nohits_seqid_file="nohits.seqid"
nohits_fasta_file="nohits.fasta"



# there are linebreaks in the OTU fasta files. Remove them
for OTU_fasta in "${parent_dir}"/*/7_OTUs.fasta; do
	awk '/^>/{print (NR==1)?$0:"\n"$0;next}{printf "%s", $0}END{print ""}' "${OTU_fasta}" > "${OTU_fasta%/*}"/OTUs.fasta	
done









for xml_file in "${parent_dir}"/*/*.xml ; do
 	N_found=$( grep -c "${string_to_match}" "${xml_file}" )
 	echo Found "${N_found}" occurrences of \""${string_to_match}"\" in "${xml_file}"
 	grep "${string_to_match}" "${xml_file}" -B 15 |  grep ''"${seq_ID_prefix}"'.*;' -o | sed 's/;//' > "${xml_file%/*}"/"${nohits_seqid_file}"
# 	echo "${xml_file%/*}"/"${nohits_seqid_file}"
	grep -f "${xml_file%/*}"/"${nohits_seqid_file}" "${xml_file%/*}"/OTUs.fasta -A 1 | sed '/--/d' > "${xml_file%/*}"/"${nohits_fasta_file}"
done


mkdir "${parent_dir}"/nohits
for tag in "${parent_dir}"/tag* ; do
# 	echo "${tag##*/}"
	cp "${tag}"/nohits.fasta "${parent_dir}"/nohits/"${tag##*/}"_nohits.fasta
done


################################################################################
# GRAVEYARD
################################################################################


# working example:
# grep 'No hits found' /Users/threeprime/Desktop/Reanalysis_20141031_1601/tag_2/8_BLASTed.xml -B 15 > temp.txt

# working version: grep 'M00.*<' temp.txt --color
# grep "${seq_ID_prefix}.*<" 

# find "${parent_dir}" -type f -name '*.xml' -exec grep -c "${string_to_match}" {} +

# for xml_file in "${xml_files[@]}"
# do
#  	N_found=$( grep -c "${string_to_match}" "${xml_file}" )
#  	echo Found "${N_found}" occurrences of \""${string_to_match}"\" in "${xml_file}"
# done
