#!/bin/bash

# Do some latter steps of the pipeline in isolation of the others

# This script should be in the same directory as the analysis_pipeline.sh file. This will find the path to that file:
SCRIPT_DIR="$(dirname "$0")"

# What is the existing directory containing the demultiplexed tag folders?

# I don't actually think we need this...
# Assign everything before the last forward slash as the path to the analysis directory, which contains
# ANALYSIS_DIR="${EXISTING_DEMULTIPLEXED_DIR}%/*"

START_TIME=$(date +%Y%m%d_%H%M)
REANALYSIS_DIR="${EXISTING_DEMULTIPLEXED_DIR%/*}"/Reanalysis_"${START_TIME}"
mkdir "${REANALYSIS_DIR}"

source "$SCRIPT_DIR/pipeline_params.sh"
cp "$SCRIPT_DIR/pipeline_params.sh" "${REANALYSIS_DIR}"/pipeline_params.txt



# Start the loop on each directory of demultiplexed tag folders.
for i in $( ls ${EXISTING_DEMULTIPLEXED_DIR} )
# ; do  echo "Processing ${EXISTING_DEMULTIPLEXED_DIR}/$i"; done

do

# CLUSTERING

mkdir "${REANALYSIS_DIR}"/${i}
CURRENT_DIR="${REANALYSIS_DIR}"/${i}

if [ ! -s "${EXISTING_DEMULTIPLEXED_DIR}/$i"/6_nosingle.fasta ]; then
	gunzip "${EXISTING_DEMULTIPLEXED_DIR}/$i"/6_nosingle.fasta.gz
fi

# CLUSTER SEQUENCES
CLUSTER_RADIUS="$(( 100 - ${CLUSTERING_PERCENT} ))"
usearch -cluster_otus "${EXISTING_DEMULTIPLEXED_DIR}/$i"/6_nosingle.fasta -otu_radius_pct "${CLUSTER_RADIUS}" -sizein -sizeout -otus "${CURRENT_DIR}"/7_OTUs.fasta -notmatched "${CURRENT_DIR}"/7_notmatched.fasta

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

done

# if [ "$PERFORM_CLEANUP" = "YES" ]; then
# 	echo "Compressing fasta and fastq files..."
# 	find "${ANALYSIS_DIR}" -type f -name '*.fasta' -exec gzip "{}" \;
# 	find "${ANALYSIS_DIR}" -type f -name '*.fastq' -exec gzip "{}" \;
# 	echo "Cleanup performed."
# else
# 	echo "Cleanup not performed."
# fi
