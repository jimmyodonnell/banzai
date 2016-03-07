#!/usr/bin/env bash

my_dir="/Users/threeprime/Desktop/20150717/libraries/quality_reports"
outfile="${my_dir}"/fastqc_counts.txt


my_files=$( find "${my_dir}" -type f -name '*.zip' )

for i in $my_files; do
	unzip $i -d "${my_dir}/temp"	
done

temp_dir=$( find "${my_dir}/temp" -type d -maxdepth 1 )

for i in $temp_dir; do
	cat "${i}"/fastqc_data.txt \
	| awk -F'\t' '/Filename|Total\ Sequences/ {print $2}' >> tmp.txt
done

sed ' {
N
s:.fastq.gz\n:,:
}' tmp.txt > "${outfile}"

rm tmp.txt

rm -rf "${my_dir}/temp"

exit



"Filename"
"Total Se"
'{for(i=1;i<=NF;i++){if(match(srch,$i)){val=$(i+1)}}}
END{print srch" "val}' input_file

Total Sequences