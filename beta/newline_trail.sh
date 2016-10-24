#!/usr/bin/env bash

# test for a trailing newline
# a trailing newline is a "return" at the end of the last line of the file

# INFILE="${1}"

if ! [ -a "${INFILE}" ]; then

  echo "Whoa there! Could not find the file you specified. Check the path."

  exit 1

fi

if [[ $(tail -c1 "${INFILE}") ]]; then

  echo "File does not have a trailing newline. I'll fix that for you."

  echo ''>>"${INFILE}"

else

  echo "File has trailing newline. Très bien!"

fi
