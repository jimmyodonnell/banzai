#!/usr/bin/env bash

INFILE="${1}"

if ! [[ -s "${INFILE}" ]]; then

  echo "Whoa there! Could not find the file you specified. Check the path."

  exit 1

fi

if [[ $( file "${INFILE}" ) == *"CRLF"* ]]; then

  echo "The file has CRLF endings. Let me fix that for you..."

  BASE="${INFILE%.*}"

  EXT="${INFILE##*.}"

  NEWLINES_FIXED="${BASE}"_fix."${EXT}"

  tr -d '\r' < "${INFILE}" > "${NEWLINES_FIXED}"

  echo "The new file is here:"

  echo "${NEWLINES_FIXED}"

else

  echo "The file passes test for CRLF. Everybody dance!"
  echo

fi
