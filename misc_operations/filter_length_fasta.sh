#!/usr/bin/env bash

# using awk

infile="${1}"
outfile=${infile%.*}

min_length="80"

awk -v min="${min_length}" '\
  BEGIN {
    RS = ">" ; ORS = ""
  } \
  length($2) >= min { \
    print ">"$0 \
  }' "${infile}" > "${outfile}"
