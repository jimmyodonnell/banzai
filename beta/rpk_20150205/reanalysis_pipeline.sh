#!/bin/bash

# Do some latter steps of the pipeline in isolation of the others

# This script should be in the same directory as the analysis_pipeline.sh file. This will find the path to that file:
SCRIPT_DIR="$(dirname "$0")"

# What is the existing directory containing the demultiplexed tag folders?

# I don't actually think we need this...
# Assign everything before the last forward slash as the path to the analysis directory, which contains
# ANALYSIS_DIR="${EXISTING_DEMULTIPLEXED_DIR}%/*"

source "$SCRIPT_DIR/pipeline_params.sh"

START_TIME=$(date +%Y%m%d_%H%M)
REANALYSIS_DIR="${EXISTING_DEMULTIPLEXED_DIR%/*}"/Reanalysis_"${START_TIME}"
mkdir "${REANALYSIS_DIR}"
cp "$SCRIPT_DIR/pipeline_params.sh" "${REANALYSIS_DIR}"/pipeline_params.txt

N_CORES=$(sysctl -n hw.ncpu)
if [ $N_CORES -gt 1 ]; then
	echo "$N_CORES cores detected."
else
	N_CORES=1
	echo "Multiple cores not detected."
fi



# Start the loop on each directory of demultiplexed tag folders.
# for i in $( ls ${EXISTING_DEMULTIPLEXED_DIR} )
# 	do
# 
# 	# CLUSTERING
# 
# 	mkdir "${REANALYSIS_DIR}"/${i}
# 	CURRENT_DIR="${REANALYSIS_DIR}"/${i}
# 
# # 	if [ ! -s "${EXISTING_DEMULTIPLEXED_DIR}/$i"/6_nosingle.fasta ]; then
# # 		gunzip "${EXISTING_DEMULTIPLEXED_DIR}/$i"/6_nosingle.fasta.gz
# # 	fi
# 
# 	# CLUSTER SEQUENCES
# 	if [ "$BLAST_WITHOUT_CLUSTERING" = "NO" ]; then
# 		CLUSTER_RADIUS="$(( 100 - ${CLUSTERING_PERCENT} ))"
# 		usearch -cluster_otus "${EXISTING_DEMULTIPLEXED_DIR}/$i"/no_duplicates.fasta -otu_radius_pct "${CLUSTER_RADIUS}" -sizein -sizeout -otus "${CURRENT_DIR}"/7_OTUs.fasta -notmatched "${CURRENT_DIR}"/7_notmatched.fasta
# 		# BLAST CLUSTERS
#		blastn -query "${CURRENT_DIR}"/7_OTUs.fasta -db "$BLAST_DB" -num_threads "$N_CORES" -perc_identity "${PERCENT_IDENTITY}" -word_size "${WORD_SIZE}" -evalue "${EVALUE}" -max_target_seqs "${MAXIMUM_MATCHES}" -outfmt 5 -out "${DEREP_INPUT%/*}"/10_BLASTed.xml

# 	else
# 		# BLAST READS
#		blastn -query "${EXISTING_DEMULTIPLEXED_DIR}/$i"/no_duplicates.fasta -db "$BLAST_DB" -num_threads "$N_CORES" -perc_identity "${PERCENT_IDENTITY}" -word_size "${WORD_SIZE}" -evalue "${EVALUE}" -max_target_seqs "${MAXIMUM_MATCHES}" -outfmt 5 -out "${DEREP_INPUT%/*}"/10_BLASTed.xml
# 	fi
# 
# done

if [ "$CONCATENATE_SAMPLES" = "YES" ]; then
	DIRECTORIES="${DEREP_INPUT%/*}"
else
	DIRECTORIES=$( find "${REANALYSIS_DIR}" -type d -d 1 )
fi


# for DIR in $DIRECTORIES; do
# 
#  		BLAST_XML=$(find "${DIR}"/10_BLASTed.xml)
#  		# Some POTENTIAL OPTIONS FOR MEGAN EXPORT:
#  		# {readname_taxonname|readname_taxonid|readname_taxonpath|readname_matches|taxonname_count|taxonpath_count|taxonid_count|taxonname_readname|taxonpath_readname|taxonid_readname}
#  		# PERFORM COMMON ANCESTOR GROUPING IN MEGAN
#  		MEGAN_COMMAND_FILE="${DIR}"/megan_commands.txt
#  		MEGAN_RMA_FILE="${DIR}"/meganfile.rma
#  		MEGAN_SHELL_SCRIPT="${DIR}"/megan_script.sh
# 
# 		#remove these files, if already present, to prevent appending
#  		rm -f ${MEGAN_COMMAND_FILE}
#  		rm -f ${MEGAN_RMA_FILE}
#  		rm -f ${MEGAN_SHELL_SCRIPT}
# 		
#  		echo "import blastfile='${BLAST_XML}' meganfile='${MEGAN_RMA_FILE}' \
#  minScore=${MINIMUM_SCORE} \
#  maxExpected=${MAX_EXPECTED} \
#  topPercent=${TOP_PERCENT} \
#  minSupportPercent=${MINIMUM_SUPPORT_PERCENT} \
#  minSupport=${MINIMUM_SUPPORT} \
#  minComplexity=${MINIMUM_COMPLEXITY} \
#  lcapercent=${LCA_PERCENT};" >> "${MEGAN_COMMAND_FILE}"
#  		echo "update;" >> "${MEGAN_COMMAND_FILE}"
#  		echo "collapse rank='$COLLAPSE_RANK1';" >> "${MEGAN_COMMAND_FILE}"
#  		echo "update;" >> "${MEGAN_COMMAND_FILE}"
#  		echo "select nodes=all;" >> "${MEGAN_COMMAND_FILE}"
#  		echo "export what=DSV format=readname_taxonname separator=comma file=${DIR}/meganout_${COLLAPSE_RANK1}.csv;" >> "${MEGAN_COMMAND_FILE}"
#  		if [ "$PERFORM_SECONDARY_MEGAN" = "YES" ]; then
#  			echo "collapse rank='$COLLAPSE_RANK2';" >> "${MEGAN_COMMAND_FILE}"
#  			echo "update;" >> "${MEGAN_COMMAND_FILE}"
#  			echo "select nodes=all;" >> "${MEGAN_COMMAND_FILE}"
#  			echo "export what=DSV format=readname_taxonname separator=comma file=${DIR}/meganout_${COLLAPSE_RANK2}.csv;" >> "${MEGAN_COMMAND_FILE}"
#  		fi
#  		echo "quit;" >> "${MEGAN_COMMAND_FILE}"
#  
#  		echo "#!/bin/bash" >> "$MEGAN_SHELL_SCRIPT"
#  #		echo "cd "${megan_exec%/*}"" >> "$MEGAN_SHELL_SCRIPT"
#  		echo "MEGAN -g -E -c ${DIR}/megan_commands.txt" >> "$MEGAN_SHELL_SCRIPT"
# 
# done
# 
# wait
# 
# for DIR in $DIRECTORIES; do
# 
# 		# Run MEGAN
# 		sh "${DIR}"/megan_script.sh
# 
# 		# Modify the MEGAN output so that it is a standard CSV file with clusterID, N_reads, and Taxon
# 		sed 's|;size=|,|' <"${DIR}"/meganout_${COLLAPSE_RANK1}.csv >"${DIR}"/meganout_${COLLAPSE_RANK1}_mod.csv
# 		sed 's|;size=|,|' <"${DIR}"/meganout_${COLLAPSE_RANK2}.csv >"${DIR}"/meganout_${COLLAPSE_RANK2}_mod.csv
# 
# 		# Run the R script, passing the current tag directory as the directory to which R will "setwd()"
# 		Rscript "$SCRIPT_DIR/megan_plotter.R" "${DIR}"
# 
# done

	
# Run the R script to consolidate samples into an OTU table; second parameter gives the relevant taxonomic level.

Rscript "$SCRIPT_DIR/consolidate_samples_020415.R" "${REANALYSIS_DIR}" "Genus"





# if [ "$PERFORM_CLEANUP" = "YES" ]; then
# 	if [ "$PIGZ_INSTALLED" = "YES" ]; then
# 		ZIPPER="pigz"
# 	else
# 		ZIPPER="gzip"
# 	fi
# 	echo "Compressing fasta and fastq files..."
# 	find "${ANALYSIS_DIR}" -type f -name '*.fasta' -exec ${ZIPPER} "{}" \;
# 	find "${ANALYSIS_DIR}" -type f -name '*.fastq' -exec ${ZIPPER} "{}" \;
# 	echo "Cleanup performed."
# else
# 	echo "Cleanup not performed."
# fi
# 
# FINISH_TIME=$(date +%Y%m%d_%H%M)
# curl http://textbelt.com/text -d number=$PHONE_NUMBER -d message="Pipeline finished! Started $START_TIME Finished $FINISH_TIME"
# 
# 
##for emailing message at completion:
LOGNAME="eDNA_Pipeline"
MAILTO="invertdna@gmail.com"
SUBJECT="Pipeline completed!"
ATTFILE=$(echo $(ls -t /Users/rpk/Desktop/*.csv | head -1))

echo -e "From: $LOGNAME\nTo: $MAILTO\nSubject: $SUBJECT\n\
Mime-Version: 1.0\nContent-Type: text/plain\n\nSample_Directory_contents:\n\n\n$(ls -lh $(echo $DIRECTORIES |cut -d' ' -f1))\n\n\nTop_Ten_Annotations\n" > /tmp/file
cat "/tmp/temp.txt" $(echo -e "\n\n") $ATTFILE>> /tmp/file
/usr/sbin/sendmail -t -oi < /tmp/file
