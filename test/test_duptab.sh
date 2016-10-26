#!/usr/bin/env bash

test_duptab () {

  file1="${1}"

  dir2="${2}"

  file2="${dir2}"/all_lib/duplicate_table.csv

  if cmp -s "$file1" "$file2"; then
    echo "The files match"
    return 0
  else
    echo "The files are different"
    return 1
  fi

}
