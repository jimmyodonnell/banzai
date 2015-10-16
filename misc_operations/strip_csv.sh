#!/usr/bin/env bash

infile="${1}"
# infile="/Users/jimmy.odonnell/Desktop/Analysis_20151013_1719/all_lib/duplicate_table.csv"
# infile="/Users/jimmy.odonnell/Desktop/toy.csv"

rownames="${infile%.*}".rownames
colnames="${infile%.*}".colnames
matrix="${infile%.*}".matrix

# awk -F',' ' { print length($0) }' "${infile}"
N_col=($( awk -F',' ' { print NF }' "${infile}" | sort -n | uniq ))
echo 'Infile has this many columns:' "${N_col[@]}"
if [[ "${#N_col[@]}" -gt 1 ]]; then
  echo 'Rows have differing number of colums. Something is wrong!'
fi

head -n 1 "${infile}" > "${colnames}"
echo 'column names written to file:' $colnames

awk -F',' ' { print $1 } ' "${infile}" > "${rownames}"
echo 'row names written to file:' $rownames

awk -F, ' NR > 1 { print $0 }' "${infile}" | cut -d',' -f 2-  > "${matrix}"
echo 'matrix values written to file:' $matrix

# awk -F',' '{print $1 > "rownames"; sub(/^[^|]+\,/,"")}1' "${infile}" > tmp && mv tmp "${matrix}"
