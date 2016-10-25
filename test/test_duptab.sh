#!/usr/bin/env bash

test_duptab () {

  dir1="${1}"

  dir2="${2}"

  file1="${dir1}"/all_lib/duplicate_table.csv

  file2="${dir2}"/all_lib/duplicate_table.csv

  if cmp -s "$file1" "$file2"; then
     echo "The files match"
  else
     echo "The files are different"
  fi

}
