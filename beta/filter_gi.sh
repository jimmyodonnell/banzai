#!/usr/bin/env bash

gi_file="/Users/threeprime/Documents/GoogleDrive/Kelly_Lab/Bioinformatics/20151113/gi_num.txt"

gi_list=($(cat $gi_file))

# path to the file containing two columns: 1 with gi and 1 with taxid
tax_dmp="/Users/threeprime/temp_big/NCBI/Taxonomy/gi_taxid_nucl.dmp"
# head -n 10000 "${tax_dmp}" > tax_dmp.txt
tax_top=$(head -n 10000 "${tax_dmp}")

outfile="somefile.txt"

# parallel --pipepart --block 100M -a tax_dmp.txt -k grep -f "${gi_file}" > "${outfile}"
# head "${tax_top}"

# fgrep -f "${gi_file}" "${tax_top}"

# echo double awk
# time awk '{print $1}' "${gi_file}" | awk '{s+=$1} END {print s}'
# echo
# echo single awk
# time awk '{s+=$1} END {print s}' "${gi_file}"
# echo
# echo cat piped to parallel
time cat "${gi_file}" | parallel --pipe awk \'{s+=\$1} END {print s}\' | awk '{s+=$1} END {print s}'
# echo
# echo no pipe
# parallel --pipepart awk \'{s+=\$1} END {print s}\'  | awk '{s+=$1} END {print s}'
exit

while IFS= read -r this_gi; do
    echo this is the line with "${this_gi}"
done < "${gi_file}"

exit

for gi in "${gi_list[@]}"; do
    echo $gi
done

exit

awk '/^'"$gi"'\t/ { print $0 }' >>


awk ' $1 == '' { print ; exit } ' $tax_dmp >>
