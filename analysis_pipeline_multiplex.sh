#!/bin/bash

# Pipeline for analysis of MULTIPLEXED Illumina data, a la Jimmy

# This command specifies the path to the directory containing the script
SCRIPT_DIR="$(dirname "$0")"

# Get the directory containing the READ1 file and assign it to variable READ_DIR.
READ_DIR="${READ1%/*}"

# Define a variable called START_TIME
START_TIME=$(date +%Y%m%d_%H%M)

# And make a directory with that timestamp
mkdir "${READ_DIR}"/Analysis_"${START_TIME}"
ANALYSIS_DIR=$"{"${READ_DIR}"/Analysis_"${START_TIME}"}"

# Read in the parameter files
source "$SCRIPT_DIR/pipeline_params.sh"
source "$SCRIPT_DIR/pear_params.sh"
CLUSTER_RADIUS="$(( 100 - ${CLUSTERING_PERCENT} ))"

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
	LINES=$(wc -l < "${ANALYSIS_DIR}"/1_merged.assembled.fastq)
	HALF_LINES=$((LINES / 2))
	head -n ${HALF_LINES} "${ANALYSIS_DIR}"/1_merged.assembled.fastq > "${ANALYSIS_DIR}"/1_merged.assembled_A.fastq
	tail -n ${HALF_LINES} "${ANALYSIS_DIR}"/1_merged.assembled.fastq > "${ANALYSIS_DIR}"/1_merged.assembled_B.fastq
	usearch -fastq_filter "${ANALYSIS_DIR}"/1_merged.assembled_A.fastq -fastq_maxee 0.5 -fastq_minlen 75 -fastaout "${ANALYSIS_DIR}"/2_filtered_A.fasta
	usearch -fastq_filter "${ANALYSIS_DIR}"/1_merged.assembled_B.fastq -fastq_maxee 0.5 -fastq_minlen 75 -fastaout "${ANALYSIS_DIR}"/2_filtered_B.fasta
	cat "${ANALYSIS_DIR}"/2_filtered_A.fasta "${ANALYSIS_DIR}"/2_filtered_B.fasta > "${ANALYSIS_DIR}"/2_filtered.fasta
else
	usearch -fastq_filter "${ANALYSIS_DIR}"/1_merged.assembled.fastq -fastq_maxee 0.5 -fastq_minlen 75 -fastaout "${ANALYSIS_DIR}"/2_filtered.fasta
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

# remove the tag from the beginning of the sequence (5' end) in it's current orientation
cutadapt -g ^NNN"${TAG}" -e 0 --discard-untrimmed "${ANALYSIS_DIR}"/3_no_homopolymers.fasta > "${CURRENT_DIR}"/5prime_tag_rm.fasta

# This sed command looks really f***ing ugly; but I'm pretty sure it works.
sed -E 's/'"${TAG}"'.{0,10}$/bananafruits/' 2_filtered.fasta > temp1.fasta

cutadapt -g ^"${PRIMER1_NON}" -g ^"${PRIMER2_NON}" -e "${PRIMER_MISMATCH_PROPORTION}" --discard-untrimmed "${CURRENT_DIR}"/5prime_tag_rm.fasta > "${CURRENT_DIR}"/5prime_primer_rm.fasta

cutadapt -a "${PRIMER1RC_NON}" -a "${PRIMER2RC_NON}" -e "${PRIMER_MISMATCH_PROPORTION}" --discard-untrimmed "${CURRENT_DIR}"/5prime_primer_rm.fasta > "${CURRENT_DIR}"/both_primer_rem.fasta




# REMOVE PRIMER SEQUENCES
cutadapt -g ^${PRIMER1_NON} -e ${PRIMER_MISMATCH_PROPORTION} --match-read-wildcards --discard-untrimmed "${ANALYSIS_DIR}"/3_no_homopolymers.fasta > "${ANALYSIS_DIR}"/PRIMER1_remtmp.fasta
cutadapt -a ${PRIMER2RC_NON} -e ${PRIMER_MISMATCH_PROPORTION} --match-read-wildcards --discard-untrimmed "${ANALYSIS_DIR}"/PRIMER1_remtmp.fasta > "${ANALYSIS_DIR}"/PRIMER1_rem.fasta
cutadapt -g ^${PRIMER2_NON} -e ${PRIMER_MISMATCH_PROPORTION} --match-read-wildcards --discard-untrimmed "${ANALYSIS_DIR}"/3_no_homopolymers.fasta > "${ANALYSIS_DIR}"/PRIMER2_remtmp.fasta
cutadapt -a ${PRIMER1RC_NON} -e ${PRIMER_MISMATCH_PROPORTION} --match-read-wildcards --discard-untrimmed "${ANALYSIS_DIR}"/PRIMER2_remtmp.fasta > "${ANALYSIS_DIR}"/PRIMER2_rem.fasta
seqtk seq -r "${ANALYSIS_DIR}"/PRIMER2_rem.fasta > "${ANALYSIS_DIR}"/PRIMER2_rc.fasta
cat "${ANALYSIS_DIR}"/PRIMER1_rem.fasta "${ANALYSIS_DIR}"/PRIMER2_rc.fasta > "${ANALYSIS_DIR}"/4_noprimers_sameorientation.fasta

# CONSOLIDATE IDENTICAL SEQUENCES. With macqiime, use: "${ANALYSIS_DIR}"/split_lib/seqs.fna
usearch -derep_fulllength "${ANALYSIS_DIR}"/4_noprimers_sameorientation.fasta -sizeout -strand both -output "${ANALYSIS_DIR}"/5_derep.fasta #-minseqlength 75

# REMOVE SINGLETONS
usearch -sortbysize "${ANALYSIS_DIR}"/5_derep.fasta -minsize 2 -sizein -sizeout -output "${ANALYSIS_DIR}"/6_nosingle.fasta

# CLUSTER SEQUENCES
usearch -cluster_otus "${ANALYSIS_DIR}"/6_nosingle.fasta -otu_radius_pct "${CLUSTER_RADIUS}" -sizein -sizeout -otus "${ANALYSIS_DIR}"/7_OTUs.fasta -notmatched "${ANALYSIS_DIR}"/7_notmatched.fasta

# BLAST CLUSTERS
blastn -query "${ANALYSIS_DIR}"/7_OTUs.fasta -db "$BLAST_DB" -perc_identity "${PERCENT_IDENTITY}" -word_size "${WORD_SIZE}" -evalue "${EVALUE}" -max_target_seqs "${MAXIMUM_MATCHES}" -outfmt 5 -out "${ANALYSIS_DIR}"/8_BLASTed.xml

# PERFORM COMMON ANCESTOR GROUPING IN MEGAN
cat > "${ANALYSIS_DIR}"/megan_commands.txt <<EOF
import blastfile='${ANALYSIS_DIR}/8_BLASTed.xml' meganfile='${ANALYSIS_DIR}/meganfile.rma';
recompute minsupport=1 mincomplexity=0;
uncollapse nodes=all;
update;
collapse rank=Family;
select nodes=all;
export what=DSV format=readname_taxonname separator=comma file=${ANALYSIS_DIR}/meganout.csv;
quit;
EOF

cat > "${ANALYSIS_DIR}"/megan_script.sh <<EOF
#!/bin/bash
cd "${megan_exec%/*}"
./"${megan_exec##*/}" -g -E -c ${ANALYSIS_DIR}/megan_commands.txt
EOF

sh "${ANALYSIS_DIR}"/megan_script.sh

# Modify the MEGAN output so that it is a standard CSV file with cluterID, N_reads, and Taxon
sed 's|;size=|,|' <"${ANALYSIS_DIR}"/meganout.csv >"${ANALYSIS_DIR}"/meganout_mod.csv

# Copy the plotting script to the current directory
cp "${SCRIPT_DIR}"/megan_plotter.R "${ANALYSIS_DIR}"/megan_plotter.R

# Add a line (before line 4) to change R's directory to the one the loop is working in (the variable ${ANALYSIS_DIR})
sed "4i\\
setwd('"${ANALYSIS_DIR}"')
" ${ANALYSIS_DIR}/megan_plotter.R > tmpfile ; mv tmpfile ${ANALYSIS_DIR}/megan_plotter.R

Rscript "${ANALYSIS_DIR}"/megan_plotter.R

if [ "$PERFORM_CLEANUP" = "YES" ]; then
	rm test1 test2 test3 test4
	echo "Cleanup performed."
else
	echo "Cleanup not performed."
fi

done
