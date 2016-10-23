#!/bin/bash

# make an otu table based on the duplicate table and usearch clustering.

MYDIR='/Users/threeprime/Desktop/Analysis_20150416_1954/all_lib'
INFILE="${MYDIR}"/'9_clusters.uc'
DUP_TABLE="${MYDIR}"/'dups.csv'
FASTA_FILE="${MYDIR}"/'no_duplicates.fasta'

# find all the sequences that are "not chimaeras" (I now can't remember what this meant...)
awk '{ if ($1 != "C") print $9 $10}' "${INFILE}" | awk -F';' '{ print $1 }' | sort -k1.5 -n | uniq | sed '/^$/d' > "${FASTA_FILE%/*}"/not_chimaeras.txt
# deciphering: here's an example line from the INFILE
# H	3	115	99.1	+	0	0	115M	DUP_8;size=5434;	DUP_4;size=14729;
# if the first field is NOT "C" (and thus either H (hit) or S (centroid)), print field 9 and 10
# e.g. DUP_8;size=5434;DUP_4;size=14729;
# then, using ';' as the field separator, print the first field (e.g. 
# e.g. DUP_8
# Sort (sort) them numerically (-n) based on the 5th character in the first field (-k1.5)
# take only the unique ones (exclude duplicates) (uniq)
# remove empty lines (sed '/^$/d')

# write a file of all of the duplicate names:
# from the duplicate table, ignore the first line (NR>1) and print just the first field (field separator defined as comma)
awk -F',' '{ if (NR > 1) print $1 }' "${DUP_TABLE}" > "${FASTA_FILE%/*}"/dup_names.txt



# 
awk '{ if ($1 == "H") print $9 $10 }' "${INFILE}" | awk -F';' '{ print $1 "," $3 }' | sort | uniq > "${FASTA_FILE%/*}"/dup_to_OTU.csv
# Code breakdown:
# If the sequence is clustered to another (a 'hit'), take fields 9 and 10 


# identify the sequences that are chimaeras (again, not sure what I meant by this)
diff "${FASTA_FILE%/*}"/dup_names.txt "${FASTA_FILE%/*}"/not_chimaeras.txt | awk '{ if ($1 == "<") print $2 }' > "${FASTA_FILE%/*}"/chimaeras.txt


awk -F'[\t;]' 'BEGIN{ print "Query,Match" } { if ($1 == "S") {print $9 "," $9} else if ($1 == "H") print $9 "," $12 }' $INFILE
# breakdown:
# When reading in the file, consider both tabs (\t) and semicolons (;) as field separators (-F)
# first print a line that reads "Query,Match"
# then read through each line of the input file, and...
# if the first field ($1) is "S", write the 9th field ($9), then a comma (","), then the 9th field again.
# if the first field ($1) is "H", write the 9th field ($9), then a comma (","), then the 12th field.