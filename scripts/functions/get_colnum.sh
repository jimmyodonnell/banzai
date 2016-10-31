#!/usr/bin/env bash

get_colnum () {
	# define a function to find column number; $1 = column name; $2 = file name
	colnum=$(awk -F',' -v COLNAME=$1 \
	  '{for (i=1;i<=NF;i++)
		    if($i == COLNAME)
			  print i;
			exit}' $2)
	if [[ "${colnum}" > 0 ]]; then
	  echo "${colnum}"
	else
		echo "ERROR! Could not find column: '""${1}""' in metadata file. Exiting."
		return 1
		exit
	fi
}
