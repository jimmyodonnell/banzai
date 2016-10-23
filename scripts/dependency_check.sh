#!/usr/bin/env bash

# check for dependencies
# usage: bash dependency_check.sh program1 program2 program3 ... programN

dependecies="$@"

echo 'Checking for dependencies:' "${dependencies[@]}"
for i in "${dependencies[@]}"; do
	if hash "${i}" 2>/dev/null; then
	# if command -v "${i}" >/dev/null 2>&1; then
		echo 'Found program' "${i}" 'in' $( which "${i}" )
	else
		echo 'ERROR: A program on which this script depends was not found:' "${i}"
		echo 'Aborting script.'
		exit
	fi
done
