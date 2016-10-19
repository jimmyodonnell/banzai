#!/usr/bin/env bash

# INFILE="${1}"

if ! [ -a "${INFILE}" ]; then

  echo "Whoa there! Could not find the file you specified. Check the path."

  exit 1

fi

if [[ $( file "${INFILE}" ) == *"CRLF"* ]]; then

  echo "The file has CRLF endings. Let me fix that for you..."

  BASE="${INFILE%.*}"

  EXT="${INFILE##*.}"

  OUTFILE="${BASE}"_fix."${EXT}"

  tr -d '\r' < "${INFILE}" > "${OUTFILE}"

  echo "The new file is here:"

  echo "${OUTFILE}"

else

  echo "File passes test for CRLF. Everybody dance!"

fi
