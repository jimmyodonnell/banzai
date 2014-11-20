#!/bin/bash

# Add the string '_LIB3' to the end of every other line, starting with line 1 (odd lines).
# awk ' NR % 2 == 1 { print $0 "_LIB3"} NR % 2 ==0 {print }' INFILE > OUTFILE

# View the last few lines of a file to make sure the renaming scheme worked the whole way through.
# tail /Users/threeprime/Documents/Data/IlluminaData/16S/run_20141113_time_series/library2/Analysis_20141113_2026/concatenated/1_demult_concat_lib2.fasta

# Replace ' 1:N:0:3' with '_'
# I originally did not write the code in the pipeline to find a '3' in that position, now it should work for any digit.
# sed s/' 1:N:0:3'/_/ /Users/threeprime/Documents/Data/IlluminaData/16S/run_20141113_time_series/library3/Analysis_20141113_2055/concatenated/1_demult_concat.fasta > /Users/threeprime/Documents/Data/IlluminaData/16S/run_20141113_time_series/library3/Analysis_20141113_2055/concatenated/1_demult_concat.fasta.tmp


# LIB1='/Users/threeprime/Documents/Data/IlluminaData/16S/run_20141113_time_series/library1/Analysis_20141113_1805/concatenated/1_demult_concat_lib1.fasta'
# LIB2='/Users/threeprime/Documents/Data/IlluminaData/16S/run_20141113_time_series/library2/Analysis_20141113_2026/concatenated/1_demult_concat_lib2.fasta'
# LIB3='/Users/threeprime/Documents/Data/IlluminaData/16S/run_20141113_time_series/library3/Analysis_20141113_2055/concatenated/1_demult_concat_lib3.fasta'

# PARENT_DIR='/Users/threeprime/Documents/Data/IlluminaData/16S/run_20141113_time_series/all_libraries'

# cat "$LIB1" "$LIB2" "$LIB3" > "$PARENT_DIR"/all_lib.fasta

# python /Applications/eDNA_analysis_pipeline/dereplicate_fasta.py "$PARENT_DIR"/all_lib.fasta

# DEREP_INPUT="$PARENT_DIR"/all_lib.fasta

# COUNT DUPLICATES PER READ, REMOVE SINGLETONS
# awk -F';' '{ if (NF > 2) print NF-1 ";" $0 }' "${DEREP_INPUT}".all | sort -nr | awk -F';' '{ print ">DUP_" NR ";" $0}' > ${DEREP_INPUT%/*}/nosingle

# NOTE LIBRARY ADDED!!!!!!!!!!!!
# PRIMER_TAGS_LIB='/Users/threeprime/Documents/Data/IlluminaData/16S/run_20141113_time_series/tags_16S_lib.txt'
TAGS=$(tr '\n' ' ' < "${PRIMER_TAGS_LIB}" )

for TAG_SEQ in $TAGS; do
	( awk 'BEGIN {print "'$TAG_SEQ'" ; FS ="_tag_'${TAG_SEQ}'" } { print NF -1 }' "${DEREP_INPUT%/*}"/nosingle > "${DEREP_INPUT%/*}"/"${TAG_SEQ}".dup ) &
done
wait

# Write a csv file of the number of occurrences of each duplicate sequence per tag.
find "${DEREP_INPUT%/*}" -type f -name '*.dup' -exec paste -d, {} \+ | awk '{ print "DUP_" NR-1 "," $0 }' > "${DEREP_INPUT%/*}"/dups.csv
rm ${DEREP_INPUT%/*}/*.dup

# Write fasta file in order to blast sequences
awk -F';' '{ print $1 ";size=" $2 ";\n" $3 }' ${DEREP_INPUT%/*}/nosingle > ${DEREP_INPUT%/*}/no_duplicates.fasta

# BLAST_INPUT="${DEREP_INPUT%/*}"/no_duplicates.fasta
# N_CORES=$(sysctl -n hw.ncpu)
# blastn -query "${BLAST_INPUT}" -db "$BLAST_DB" -num_threads "$N_CORES" -perc_identity "${PERCENT_IDENTITY}" -word_size "${WORD_SIZE}" -evalue "${EVALUE}" -max_target_seqs "${MAXIMUM_MATCHES}" -outfmt 5 -out "${DEREP_INPUT%/*}"/10_BLASTed.xml
# BLAST_XML="${DEREP_INPUT%/*}"/10_BLASTed.xml
