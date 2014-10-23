#!/bin/bash

# Pipeline for analysis of MULTIPLEXED Illumina data, a la Jimmy

# An attempt to cause the script to exit if any of the commands returns a non-zero status (i.e. FAILS).
# set -e

# This command specifies the path to the directory containing the script
SCRIPT_DIR="$(dirname "$0")"

# Read in the parameter files
source "$SCRIPT_DIR/pipeline_params.sh"
source "$SCRIPT_DIR/pear_params.sh"

# Get the directory containing the READ1 file and assign it to variable READ_DIR.
READ_DIR="${READ1%/*}"

# Define a variable called START_TIME
START_TIME=$(date +%Y%m%d_%H%M)

# And make a directory with that timestamp
mkdir "${READ_DIR}"/Analysis_"${START_TIME}"
ANALYSIS_DIR="${READ_DIR}"/Analysis_"${START_TIME}"

# Copy these files into that directory as a verifiable log you can refer back to.
cp "${SCRIPT_DIR}"/analysis_pipeline_multiplex.sh "${ANALYSIS_DIR}"/analysis_pipeline_used.txt
cp "${SCRIPT_DIR}"/pipeline_params.sh "${ANALYSIS_DIR}"/pipeline_parameters.txt
cp "${SCRIPT_DIR}"/pear_params.sh "${ANALYSIS_DIR}"/pear_parameters.txt

# Read in primers and their reverse complements.
PRIMER1=$( awk 'NR==2' "${PRIMER_FILE}" )
PRIMER2=$( awk 'NR==4' "${PRIMER_FILE}" )
PRIMER1RC=$( seqtk seq -r "${PRIMER_FILE}" | awk 'NR==2' )
PRIMER2RC=$( seqtk seq -r "${PRIMER_FILE}" | awk 'NR==4' )

# Take ambiguities out of primers. The sed command says "turn any character that's not A, T, C, or G, and replace it with N.
PRIMER1_NON=$( echo $PRIMER1 | sed "s/[^ATCG]/N/g" )
PRIMER2_NON=$( echo $PRIMER2 | sed "s/[^ATCG]/N/g" )
PRIMER1RC_NON=$( echo $PRIMER1RC | sed "s/[^ATCG]/N/g" )
PRIMER2RC_NON=$( echo $PRIMER2RC | sed "s/[^ATCG]/N/g" )

#MERGE PAIRED-END READS (PEAR)
pear -f "${READ1}" -r "${READ2}" -o "${ANALYSIS_DIR}"/1_merged -v $MINOVERLAP -m $ASSMAX -n $ASSMIN -t $TRIMMIN -q $QT -u $UNCALLEDMAX -g $TEST -p $PVALUE -s $SCORING -j $THREADS

# FILTER READS (This is the last step that uses quality scores, so convert to fasta)
# The 32bit version of usearch will not accept an input file greater than 4GB. The 64bit usearch is $900. Thus, for now:
INFILE_SIZE=$(stat "${ANALYSIS_DIR}"/1_merged.assembled.fastq | awk '{ print $8 }')
if [ ${INFILE_SIZE} -gt 4000000000 ]; then
# Must first check the number of reads. If odd, file must be split so as not to split the middle read's sequence from its quality score.
	LINES_MERGED=$(wc -l < "${ANALYSIS_DIR}"/1_merged.assembled.fastq)
	READS_MERGED=$(( LINES_MERGED / 4 ))
	HALF_LINES=$((LINES_MERGED / 2))
	if [ $((READS_MERGED%2)) -eq 0 ]; then
		head -n ${HALF_LINES} "${ANALYSIS_DIR}"/1_merged.assembled.fastq > "${ANALYSIS_DIR}"/1_merged.assembled_A.fastq
		tail -n ${HALF_LINES} "${ANALYSIS_DIR}"/1_merged.assembled.fastq > "${ANALYSIS_DIR}"/1_merged.assembled_B.fastq
	else
		head -n $(( HALF_LINES + 2 )) "${ANALYSIS_DIR}"/1_merged.assembled.fastq > "${ANALYSIS_DIR}"/1_merged.assembled_A.fastq
		tail -n $(( HALF_LINES - 2 )) "${ANALYSIS_DIR}"/1_merged.assembled.fastq > "${ANALYSIS_DIR}"/1_merged.assembled_B.fastq
	fi
	usearch -fastq_filter "${ANALYSIS_DIR}"/1_merged.assembled_A.fastq -fastq_maxee 0.5 -fastq_minlen $ASSMIN -fastaout "${ANALYSIS_DIR}"/2_filtered_A.fasta
	usearch -fastq_filter "${ANALYSIS_DIR}"/1_merged.assembled_B.fastq -fastq_maxee 0.5 -fastq_minlen $ASSMIN -fastaout "${ANALYSIS_DIR}"/2_filtered_B.fasta
	cat "${ANALYSIS_DIR}"/2_filtered_A.fasta "${ANALYSIS_DIR}"/2_filtered_B.fasta > "${ANALYSIS_DIR}"/2_filtered.fasta
else
	usearch -fastq_filter "${ANALYSIS_DIR}"/1_merged.assembled.fastq -fastq_maxee 0.5 -fastq_minlen $ASSMIN -fastaout "${ANALYSIS_DIR}"/2_filtered.fasta
fi

# REMOVE SEQUENCES CONTAINING HOMOPOLYMERS
grep -E -i "(A|T|C|G)\1{$HOMOPOLYMER_MAX,}" "${ANALYSIS_DIR}"/2_filtered.fasta -B 1 -n | cut -f1 -d: | cut -f1 -d- | sed '/^$/d' > "${ANALYSIS_DIR}"/homopolymer_line_numbers.txt
if [ -s "${ANALYSIS_DIR}"/homopolymer_line_numbers.txt ]; then
	awk 'NR==FNR{l[$0];next;} !(FNR in l)' "${ANALYSIS_DIR}"/homopolymer_line_numbers.txt "${ANALYSIS_DIR}"/2_filtered.fasta > "${ANALYSIS_DIR}"/3_no_homopolymers.fasta
	awk 'NR==FNR{l[$0];next;} (FNR in l)' "${ANALYSIS_DIR}"/homopolymer_line_numbers.txt "${ANALYSIS_DIR}"/2_filtered.fasta > "${ANALYSIS_DIR}"/homopolymeric_reads.fasta
else
	cp "${ANALYSIS_DIR}"/2_filtered.fasta "${ANALYSIS_DIR}"/3_no_homopolymers.fasta
fi

# DEMULTIPLEXING STARTS HERE
# make a directory to put all the demultiplexed files in
mkdir "${ANALYSIS_DIR}"/demultiplexed

N_TAGS=$( wc -l < "${PRIMER_TAGS}" )

# Write a file of sequence names to make a tag fasta file (necessary for reverse complementing)
for i in `seq ${N_TAGS}`; do echo \>tag"$i"; done > "${ANALYSIS_DIR}"/tag_names.txt
# Alternately paste those names and the sequences to make a tag fasta file.
paste -d"\n" "${ANALYSIS_DIR}"/tag_names.txt "${PRIMER_TAGS}" > "${ANALYSIS_DIR}"/tags.fasta
# Reverse complement the tags
seqtk seq -r "${ANALYSIS_DIR}"/tags.fasta > "${ANALYSIS_DIR}"/tags_RC.fasta

# Start the loop, do one loop for each of the number of lines in the tag file.
for (( i=1; i<=${N_TAGS}; i++ ));

do

# Assign the current tag to to variable TAG, and its reverse complement to TAG_RC
TAG=$( sed -n $((i * 2))p "${ANALYSIS_DIR}"/tags.fasta )
TAG_RC=$( sed -n $((i * 2))p "${ANALYSIS_DIR}"/tags_RC.fasta )

# Create a directory for the tag
mkdir "${ANALYSIS_DIR}"/demultiplexed/tag_"${i}"

# Make a variable (CURRENT_DIR) with the current tag's directory for ease of reading and writing.
CURRENT_DIR="${ANALYSIS_DIR}"/demultiplexed/tag_"${i}"

# REMOVE TAG SEQUENCES
# remove the tag from the beginning of the sequence (5' end) in it's current orientation
cutadapt -g ^NNN"${TAG}" -e 0 --discard-untrimmed "${ANALYSIS_DIR}"/3_no_homopolymers.fasta > "${CURRENT_DIR}"/5prime_tag_rm.fasta
# Need to first grep lines containing pattern first, THEN following sed command with remove them
grep -E "${TAG_RC}.{0,9}$" -B 1 "${CURRENT_DIR}"/5prime_tag_rm.fasta | grep -v -- "^--$"  > "${CURRENT_DIR}"/3_prime_tagged.fasta
# Turn the following line on to write chimaeras to a file.
# grep -E -v "${TAG_RC}.{0,9}$" "${CURRENT_DIR}"/5prime_tag_rm.fasta > "${CURRENT_DIR}"/chimaeras.fasta
# This sed command looks really f***ing ugly; but I'm pretty sure it works.
sed -E 's/'"${TAG_RC}"'.{0,9}$//' "${CURRENT_DIR}"/3_prime_tagged.fasta > "${CURRENT_DIR}"/3prime_tag_rm.fasta

# REMOVE PRIMER SEQUENCES
# Remove PRIMER1 and PRIMER2 from the beginning of the reads.
cutadapt -g ^"${PRIMER1_NON}" -g ^"${PRIMER2_NON}" -e "${PRIMER_MISMATCH_PROPORTION}" --discard-untrimmed "${CURRENT_DIR}"/3prime_tag_rm.fasta > "${CURRENT_DIR}"/5prime_primer_rm.fasta
# Remove the reverse complement of PRIMER1 and PRIMER2 from the end of the reads. NOTE cutadapt1.7
cutadapt -a "${PRIMER1RC_NON}" -a "${PRIMER2RC_NON}" -e "${PRIMER_MISMATCH_PROPORTION}" --discard-untrimmed "${CURRENT_DIR}"/5prime_primer_rm.fasta > "${CURRENT_DIR}"/both_primer_rem.fasta

# CONSOLIDATE IDENTICAL SEQUENCES.
usearch -derep_fulllength "${CURRENT_DIR}"/both_primer_rem.fasta -sizeout -strand both -output "${CURRENT_DIR}"/5_derep.fasta

# REMOVE SINGLETONS
usearch -sortbysize "${CURRENT_DIR}"/5_derep.fasta -minsize 2 -sizein -sizeout -output "${CURRENT_DIR}"/6_nosingle.fasta

# CLUSTER SEQUENCES
CLUSTER_RADIUS="$(( 100 - ${CLUSTERING_PERCENT} ))"
usearch -cluster_otus "${CURRENT_DIR}"/6_nosingle.fasta -otu_radius_pct "${CLUSTER_RADIUS}" -sizein -sizeout -otus "${CURRENT_DIR}"/7_OTUs.fasta -notmatched "${CURRENT_DIR}"/7_notmatched.fasta

# BLAST CLUSTERS
blastn -query "${CURRENT_DIR}"/7_OTUs.fasta -db "$BLAST_DB" -perc_identity "${PERCENT_IDENTITY}" -word_size "${WORD_SIZE}" -evalue "${EVALUE}" -max_target_seqs "${MAXIMUM_MATCHES}" -outfmt 5 -out "${CURRENT_DIR}"/8_BLASTed.xml

# PERFORM COMMON ANCESTOR GROUPING IN MEGAN
cat > "${CURRENT_DIR}"/megan_commands.txt <<EOF
import blastfile='${CURRENT_DIR}/8_BLASTed.xml' meganfile='${CURRENT_DIR}/meganfile.rma';
recompute minsupport=1 mincomplexity=0;
uncollapse nodes=all;
update;
collapse rank='$COLLAPSE_RANK';
select nodes=all;
export what=DSV format=readname_taxonname separator=comma file=${CURRENT_DIR}/meganout.csv;
quit;
EOF

cat > "${CURRENT_DIR}"/megan_script.sh <<EOF
#!/bin/bash
cd "${megan_exec%/*}"
./"${megan_exec##*/}" -g -E -c ${CURRENT_DIR}/megan_commands.txt
EOF

sh "${CURRENT_DIR}"/megan_script.sh

# Modify the MEGAN output so that it is a standard CSV file with cluterID, N_reads, and Taxon
sed 's|;size=|,|' <"${CURRENT_DIR}"/meganout.csv >"${CURRENT_DIR}"/meganout_mod.csv

# Copy the plotting script to the current directory
cp "${SCRIPT_DIR}"/megan_plotter.R "${CURRENT_DIR}"/megan_plotter.R

# PLOTTING ANNOTATIONS
# Add a line (before line 4) to change R's directory to the one the loop is working in (the variable ${CURRENT_DIR}), and copy to the current directory
sed "4i\\
setwd('"${CURRENT_DIR}"')
" ${SCRIPT_DIR}/megan_plotter.R > ${CURRENT_DIR}/megan_plotter.R

Rscript "${CURRENT_DIR}"/megan_plotter.R

if [ "$PERFORM_CLEANUP" = "YES" ]; then
	rm test1 test2 test3 test4
	echo "Cleanup performed."
else
	echo "Cleanup not performed."
fi

done
