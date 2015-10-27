#!/usr/bin/env bash

# usage:
bash characters_adjacent.sh myfile.txt 

infile="${1}"
pattern="_lib_"
char_following=1


# Don't edit this
pattern_length=${#pattern}
awk_rlength=$(( char_following - pattern_length ))

awk '/^>/ {match($0, /'"${pattern}"'/); print substr($0, RSTART + RLENGTH, '"${char_following}"');}' "${infile}"

# '"${pattern_length}"'

exit



# Substring works like this:
#substr(string, start, length)
# (prints to the end of the line if no "length" argument given)

mystring=abcdefghijklmnopqrstuvwxyz

# print the string and the rest of the line
echo $mystring$mystring | awk '{match($0, /def/); print substr($0, RSTART) }'


# offset the starting of the printing using RLENGTH:
echo $mystring$mystring | awk '{match($0, /def/); print substr($0, RSTART + 2) }'
# fghijklmnopqrstuvwxyzabcdefghijklmnopqrstuvwxyz

# print the entire string following the line:
echo $mystring$mystring | awk '{match($0, /def/); print substr($0, RSTART + RLENGTH) }'
# ghijklmnopqrstuvwxyzabcdefghijklmnopqrstuvwxyz


echo $mystring$mystring | awk '{match($0, /def/); print substr($0, RSTART + RLENGTH, 6) }'


# this is the demultiplexing part from the pipeline:
awk 'gsub(/^.{3}'"$TAG_SEQ"'/,"") {if (a && a !~ /^.{3}'"$TAG_SEQ"'/) print a; print} {a=$0}'