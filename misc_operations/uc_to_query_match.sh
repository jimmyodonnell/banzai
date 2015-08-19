#!/usr/bin/env bash

# convert usearch/uparse-like output, e.g.:

# to a file with only two columns, Query and Match, e.g.:
# DUP_278557,DUP_32996

# if the input file looks like this:
# H	0	169	99.4	+	0	0	169M	DUP_278557;size=4;	DUP_32996;size=48;
# set uc_type="7"

# if the uc file looks like this:
# DUP_5;size=1317763;	otu	83.2	*	OTU_1	OTU_5
# set uc_type="7"

uc_type="8"

infile="${1}"

outfile="${infile%/*}"/dups_to_otus.csv

echo "${outfile}"


if [ "$uc_type" = "7" ]

then

	awk -F'[\t;]' 'BEGIN{ print "Query,Match" } { if ($1 == "S") {print $9 "," $9} else if ($1 == "H") { print $9 "," $12 } }' "${infile}" > "${outfile}"

elif [ "$uc_type" = "8" ]
		
then
	
	awk -F'[\t;]' 'BEGIN{ print "Query,Match" } { if ($4 == "otu") {print $1 "," $1} else if ($4 == "match") { print $1 "," $7 } else if ($4 == "chimera") { print $1 "," "chimera"} }' "${infile}" > "${outfile}"
	
else
	
	echo "incorrect uc_type assignment"	
	
fi




