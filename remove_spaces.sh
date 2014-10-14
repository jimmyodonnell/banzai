#!/bin/bash

# put this file into the directory of interest.
# upon execution it will replace the spaces (' ') in directory names with underscores ('_')

mydir="$(dirname "$0")"
par_dir="${mydir%/*}"
find $par_dir -name '* *' -type d | while read file;
do
target=`echo "$file" | sed 's/ /_/g'`;
echo "Renaming '$file' to '$target'";
#mv "$file" "$target";
done;
