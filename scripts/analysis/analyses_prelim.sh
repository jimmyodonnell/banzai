#!/usr/bin/env bash

################################################################################
# PRELIMINARY ANALYSES
################################################################################
# Once you have a final CSV file of the number of occurences of each OTU in each sample, run some preliminary analyses in R
# TODO rename preliminary to OTU analyses; move analysis script to OTU analysis directory
OUTPUT_PDF="${OUTPUT_DIR}"/analysis_results_"${START_TIME}".pdf

echo $(date +%Y-%m-%d\ %H:%M) "passing args to R for preliminary analysis..."
Rscript "$SCRIPT_DIR/scripts/analysis/analyses_prelim.R" "${OUTPUT_PDF}" "${OTU_table}" "${SEQUENCING_METADATA}" "${LIBRARY_COLUMN_NAME}" "${SECONDARY_INDEX_COLUMN_NAME}" "${ColumnName_SampleName}" "${ColumnName_SampleType}"
echo

# EMPTY PDFs are 3829 bytes
minimumsize=4000
size_PDF=$(wc -c <"${OUTPUT_PDF}")
if [ "${size_PDF}" -lt "${minimumsize}" ]; then
    echo 'There was a problem generating the PDF.'
else
	REMOTE_PDF="${OUTPUT_PDF_DIR}"/analysis_results_"${START_TIME}".pdf
	cp "${OUTPUT_PDF}" "${REMOTE_PDF}"
fi
