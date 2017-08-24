#!/usr/bin/env bash

# Author: Jimmy O'Donnell <jodonnellbio@gmail.com>

################################################################################

target_dir="${1}"

if [[ ! -d "${target_dir}" ]]; then
    echo Could not find directory "${target_dir}". Exiting...
    exit 1
fi

outfile="${2}"

if [[ -f "${outfile}" ]]; then
    echo "${outfile}" exists. Overwriting...
fi

egrep '(library\(|require\()' -r "${target_dir}" \
    --include \*.r --include \*.R  |\
    awk -F'[()]' '{ print $2 }' |\
    sort |\
    uniq > "${outfile}"

echo "Required R packages:"
cat "${outfile}"

echo
echo "Written to ""${outfile}"
