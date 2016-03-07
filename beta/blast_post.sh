#!/usr/bin/env bash


blastout="/Users/jimmy.odonnell/Desktop/Analysis_20151019_1918/all_lib/OTUs_swarm/blast_20151103_1808/blasted_20151103_1808_i99.txt"


# head "${blastout}"

cat "${blastout}" | awk '{ print $4 }' | sort -n | uniq -c
exit




awk '{ print $1 }'


# print sequence length
infile="/Users/jimmy.odonnell/Desktop/Analysis_20151019_1918/all_lib/OTUs_swarm/OTUs_nochime_filt.fasta"

awk ' !/^>/ { print length } '




#!/usr/bin/env bash

cd ~/Desktop/Analysis_20151019_1918/all_lib/OTUs_swarm/blast_20151110_1606

infile="blasted_20151110_1606_e1.485037e-40.txt"

infile_seqids="blasted_20151110_1606_e1.485037e-40_nohits.names"

hits="blasted_20151110_1606_e1.031588e-36.txt"

awk '{ print $1 }' "${hits}" | uniq | wc -l
wc -l "${infile_seqids}"

#| wc -l

          # awk '{ print $1 }' "${hits}" | uniq | comm -31 - "${infile_seqids}" | wc -l


# comm -1 <(sort $infile_seqids) <(awk '{ print $1 }' "${hits}" | uniq | sort)

# grep -f <(awk '{ print $1 }' "${hits}" | uniq | sort) <(sort $infile_seqids) -c
# | comm -31 - "${infile_seqids}"

exit
