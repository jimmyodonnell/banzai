#!/usr/bin/env bash

################################################################################
## MEGAN ##
PERFORM_MEGAN="NO"

# For more information, see the manual provided with the software
# Specify the path to the MEGAN executable file you want to use.
# Note that in recent versions an executable was not provided; in that case, you need to reference like so: '/Applications/MEGAN/MEGAN.app/Contents/MacOS/JavaApplicationStub'
megan_exec='/Applications/MEGAN/MEGAN.app/Contents/MacOS/JavaApplicationStub'

# What is the lowest taxonomic rank at which MEGAN should group OTUs?
COLLAPSE_RANK1="Family"
MINIMUM_SUPPORT="1"
MINIMUM_COMPLEXITY="0"
TOP_PERCENT="3"
MINIMUM_SUPPORT_PERCENT="0"
MINIMUM_SCORE="140"
LCA_PERCENT="70"
MAX_EXPECTED="1e-25"

# Do you want to perform a secondary MEGAN analysis, collapsing at a different taxonomic level?
PERFORM_SECONDARY_MEGAN="YES"
COLLAPSE_RANK2="Genus"

if [ "$PERFORM_MEGAN" = "YES" ]; then

	DIRECTORIES="${DEREP_INPUT%/*}"

	for DIR in "$DIRECTORIES"; do

		# Some POTENTIAL OPTIONS FOR MEGAN EXPORT:
		# {readname_taxonname|readname_taxonid|readname_taxonpath|readname_matches|taxonname_count|taxonpath_count|taxonid_count|taxonname_readname|taxonpath_readname|taxonid_readname}
		# PERFORM COMMON ANCESTOR GROUPING IN MEGAN

			# check for blast output
		if [[ -s "${blast_output}"  ]]; then

			echo $(date +%Y-%m-%d\ %H:%M) 'BLAST output found; proceeding to MEGAN.'
			echo
			# Specify paths to megan-related files
			BLAST_XML="${DIR}"/10_BLASTed.xml
			MEGAN_COMMAND_FILE="${DIR}"/megan_commands.txt
			MEGAN_RMA_FILE="${DIR}"/meganfile.rma
			MEGAN_SHELL_SCRIPT="${DIR}"/megan_script.sh

			echo "import blastfile='${BLAST_XML}' meganFile='${MEGAN_RMA_FILE}' \
	minScore=${MINIMUM_SCORE} \
	maxExpected=${MAX_EXPECTED} \
	topPercent=${TOP_PERCENT} \
	minSupportPercent=${MINIMUM_SUPPORT_PERCENT} \
	minSupport=${MINIMUM_SUPPORT} \
	minComplexity=${MINIMUM_COMPLEXITY} \
	lcapercent=${LCA_PERCENT};" > "${MEGAN_COMMAND_FILE}"
			echo "update;" >> "${MEGAN_COMMAND_FILE}"
			echo "collapse rank='$COLLAPSE_RANK1';" >> "${MEGAN_COMMAND_FILE}"
			echo "update;" >> "${MEGAN_COMMAND_FILE}"
			echo "select nodes=all;" >> "${MEGAN_COMMAND_FILE}"
			echo "export what=DSV format=readname_taxonname separator=comma file=${DIR}/meganout_${COLLAPSE_RANK1}.csv;" >> "${MEGAN_COMMAND_FILE}"
			if [ "$PERFORM_SECONDARY_MEGAN" = "YES" ]; then
				echo "collapse rank='$COLLAPSE_RANK2';" >> "${MEGAN_COMMAND_FILE}"
				echo "update;" >> "${MEGAN_COMMAND_FILE}"
				echo "select nodes=all;" >> "${MEGAN_COMMAND_FILE}"
				echo "export what=DSV format=readname_taxonname separator=comma file=${DIR}/meganout_${COLLAPSE_RANK2}.csv;" >> "${MEGAN_COMMAND_FILE}"
			fi
			echo "quit;" >> "${MEGAN_COMMAND_FILE}"

			echo "#!/bin/bash" > "$MEGAN_SHELL_SCRIPT"
			echo "cd "${megan_exec%/*}"" >> "$MEGAN_SHELL_SCRIPT"
			echo "./"${megan_exec##*/}" -g -E -c ${DIR}/megan_commands.txt" >> "$MEGAN_SHELL_SCRIPT"

			# Run MEGAN
			sh "${DIR}"/megan_script.sh

			# Modify the MEGAN output so that it is a standard CSV file with clusterID, N_reads, and Taxon
			sed 's|;size=|,|' <"${DIR}"/meganout_${COLLAPSE_RANK1}.csv >"${DIR}"/meganout_${COLLAPSE_RANK1}_mod.csv
			sed 's|;size=|,|' <"${DIR}"/meganout_${COLLAPSE_RANK2}.csv >"${DIR}"/meganout_${COLLAPSE_RANK2}_mod.csv

			# Run the R script, passing the current tag directory as the directory to which R will "setwd()"
			Rscript "$SCRIPT_DIR/scripts/analysis/megan_plotter.R" "${DIR}"

		else
			echo
			echo 'BLAST failed: the output file is empty or absent.'
			echo 'File should be:' "${blast_output}"
			echo
		fi

	done

else
	echo 'MEGAN not performed.'
	echo
fi
