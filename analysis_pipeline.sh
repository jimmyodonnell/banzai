#!/bin/bash

# Pipeline for analysis of Illumina data, a la Jimmy

# This command specifies the path to the directory containing the script
SCRIPT_DIR="$(dirname "$0")"

# This command specifies the path to the directory containing the above directory, which should contain all the folders you'd like analyzed.
PAR_DIR="${SCRIPT_DIR%/*}"

# Read in the parameter files
source "$SCRIPT_DIR/pipeline_params.sh"
source "$SCRIPT_DIR/pear_params.sh"
CLUSTER_RADIUS="$(( 100 - ${CLUSTERING_PERCENT} ))"

# Read in primers and their reverse complements.
PRIMER1=$( awk 'NR==2' "$SCRIPT_DIR/primers.fasta" )
PRIMER2=$( awk 'NR==4' "$SCRIPT_DIR/primers.fasta" )
PRIMER1RC=$( seqtk seq -r "$SCRIPT_DIR/primers.fasta" | awk 'NR==2' )
PRIMER2RC=$( seqtk seq -r "$SCRIPT_DIR/primers.fasta" | awk 'NR==4' )

# Take ambiguities out of primers. The sed command says "turn any character that's not A, T, C, or G, and replace it with N.
PRIMER1_NON=$( echo $PRIMER1 | sed "s/[^ATCG]/N/g" )
PRIMER2_NON=$( echo $PRIMER2 | sed "s/[^ATCG]/N/g" )
PRIMER1RC_NON=$( echo $PRIMER1RC | sed "s/[^ATCG]/N/g" )
PRIMER2RC_NON=$( echo $PRIMER2RC | sed "s/[^ATCG]/N/g" )

# Define the list of folders to analyze by finding folders that contain fastq files in the same parent directory as this analysis pipeline folder.
FOLDERLIST=( $(find "$PAR_DIR" -name '*.fastq' -print0 | xargs -0 -n1 dirname | sort --unique) )

# Define the number of directories to analyze
N_DIR=${#FOLDERLIST[@]}

# Initiate pipeline loop!
for (( i=0; i<${N_DIR}; i++ ));

do

# For readability, instead of referring to the current directory the script is analyzing as FOLDERLIST[$i]:
# Assign the current directory to variable CURRENT_DIR. (switch back by find(CURRENT_DIR, replace(FOLDERLIST[$i]))
CURRENT_DIR=${FOLDERLIST[$i]}

# Identify the forward and reverse fastq files.
READ1=$(find "${CURRENT_DIR}" -name '*R1*fastq')
READ2=$(find "${CURRENT_DIR}" -name '*R2*fastq')

#MERGE PAIRED-END READS (PEAR) (consider piping pear output by adding > pear_log.txt to end of line)
pear -f "${READ1}" -r "${READ2}" -o "${CURRENT_DIR}"/1_merged -v $MINOVERLAP -m $ASSMAX -n $ASSMIN -t $TRIMMIN -q $QT -u $UNCALLEDMAX -g $TEST -p $PVALUE -s $SCORING -j $THREADS 

# FILTER READS (This is the last step that uses quality scores, so convert to fasta)
usearch -fastq_filter "${CURRENT_DIR}"/1_merged.assembled.fastq -fastq_maxee 0.5 -fastq_minlen 75 -fastaout "${CURRENT_DIR}"/2_filtered.fasta

# REMOVE SEQUENCES CONTAINING HOMOPOLYMERS
grep -E -i "(A|T|C|G)\1{$HOMOPOLYMER_MAX,}" "${CURRENT_DIR}"/2_filtered.fasta -B 1 -n | cut -f1 -d: | cut -f1 -d- | sed '/^$/d' > "${CURRENT_DIR}"/homopolymer_line_numbers.txt
if [ -s "${CURRENT_DIR}"/homopolymer_line_numbers.txt ]; then
	awk 'NR==FNR{l[$0];next;} !(FNR in l)' "${CURRENT_DIR}"/homopolymer_line_numbers.txt "${CURRENT_DIR}"/2_filtered.fasta > "${CURRENT_DIR}"/3_no_homopolymers.fasta
	awk 'NR==FNR{l[$0];next;} (FNR in l)' "${CURRENT_DIR}"/homopolymer_line_numbers.txt "${CURRENT_DIR}"/2_filtered.fasta > "${CURRENT_DIR}"/homopolymeric_reads.fasta
else
	cp "${CURRENT_DIR}"/2_filtered.fasta "${CURRENT_DIR}"/3_no_homopolymers.fasta
fi
# rm "${CURRENT_DIR}"/homopolymer_line_numbers.txt

# REMOVE PRIMER SEQUENCES
cutadapt -g ^${PRIMER1_NON} -e ${PRIMER_MISMATCH_PROPORTION} --match-read-wildcards --discard-untrimmed "${CURRENT_DIR}"/3_no_homopolymers.fasta > "${CURRENT_DIR}"/PRIMER1_remtmp.fasta
cutadapt -a ${PRIMER2RC_NON} -e ${PRIMER_MISMATCH_PROPORTION} --match-read-wildcards --discard-untrimmed "${CURRENT_DIR}"/PRIMER1_remtmp.fasta > "${CURRENT_DIR}"/PRIMER1_rem.fasta
cutadapt -g ^${PRIMER2_NON} -e ${PRIMER_MISMATCH_PROPORTION} --match-read-wildcards --discard-untrimmed "${CURRENT_DIR}"/3_no_homopolymers.fasta > "${CURRENT_DIR}"/PRIMER2_remtmp.fasta
cutadapt -a ${PRIMER1RC_NON} -e ${PRIMER_MISMATCH_PROPORTION} --match-read-wildcards --discard-untrimmed "${CURRENT_DIR}"/PRIMER2_remtmp.fasta > "${CURRENT_DIR}"/PRIMER2_rem.fasta
seqtk seq -r "${CURRENT_DIR}"/PRIMER2_rem.fasta > "${CURRENT_DIR}"/PRIMER2_rc.fasta
cat "${CURRENT_DIR}"/PRIMER1_rem.fasta "${CURRENT_DIR}"/PRIMER2_rc.fasta > "${CURRENT_DIR}"/4_noprimers_sameorientation.fasta

# CONSOLIDATE IDENTICAL SEQUENCES. With macqiime, use: "${CURRENT_DIR}"/split_lib/seqs.fna
usearch -derep_fulllength "${CURRENT_DIR}"/4_noprimers_sameorientation.fasta -sizeout -strand both -output "${CURRENT_DIR}"/5_derep.fasta #-minseqlength 75 

# REMOVE SINGLETONS
usearch -sortbysize "${CURRENT_DIR}"/5_derep.fasta -minsize 2 -sizein -sizeout -output "${CURRENT_DIR}"/6_nosingle.fasta

# CLUSTER SEQUENCES
usearch -cluster_otus "${CURRENT_DIR}"/6_nosingle.fasta -otu_radius_pct "${CLUSTER_RADIUS}" -sizein -sizeout -otus "${CURRENT_DIR}"/7_OTUs.fasta -notmatched "${CURRENT_DIR}"/7_notmatched.fasta

# BLAST CLUSTERS
blastn -query "${CURRENT_DIR}"/7_OTUs.fasta -db "$BLAST_DB"/*.fasta -perc_identity "${PERCENT_IDENTITY}" -word_size "${WORD_SIZE}" -max_target_seqs "${MAXIMUM_MATCHES}" -outfmt 5 -out "${CURRENT_DIR}"/8_BLASTed.xml
#/16Smetazoa.fasta
# usearch -usearch_global "${CURRENT_DIR}"/5_OTUs.fasta -db "$BLAST_DB" -strand plus -id 0.97 -uc "${CURRENT_DIR}"/readmap.uc

#python "${SCRIPT_DIR}"/uc2otutab.py "${CURRENT_DIR}"/readmap.uc > "${CURRENT_DIR}"/tabfile

# PERFORM COMMON ANCESTOR GROUPING IN MEGAN
cat > "${CURRENT_DIR}"/megan_commands.txt <<EOF
import blastfile='${CURRENT_DIR}/8_BLASTed.xml' meganfile='${CURRENT_DIR}/meganfile.rma';
recompute minsupport=1 mincomplexity=0;
uncollapse nodes=all;
update;
collapse rank=Family;
select nodes=all;
export what=DSV format=readname_taxonname separator=comma file=${CURRENT_DIR}/meganout.csv;
quit;
EOF

cat > "${CURRENT_DIR}"/megan_script.sh <<EOF
#!/bin/bash
cd ${megan_dir}
./MEGAN -g -E -c ${CURRENT_DIR}/megan_commands.txt
EOF

sh "${CURRENT_DIR}"/megan_script.sh

# Modify the MEGAN output so that it is a standard CSV file with cluterID, N_reads, and Taxon
sed 's|;size=|,|' <"${CURRENT_DIR}"/meganout.csv >"${CURRENT_DIR}"/meganout_mod.csv

# Copy the plotting script to the current directory
cp "${SCRIPT_DIR}"/megan_plotter.R "${CURRENT_DIR}"/megan_plotter.R

# Add a line (before line 4) to change R's directory to the one the loop is working in (the variable ${CURRENT_DIR})
sed "4i\\
setwd('"${CURRENT_DIR}"')
" ${CURRENT_DIR}/megan_plotter.R > tmpfile ; mv tmpfile ${CURRENT_DIR}/megan_plotter.R

Rscript "${CURRENT_DIR}"/megan_plotter.R

if [ PERFORM_CLEANUP=="YES" ]; then
	rm test1 test2 test3 test4
	echo "Cleanup performed."
fi

done