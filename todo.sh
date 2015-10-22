#!/usr/bin/env bash

# TODO add flexibility to stripchart "Reads by sample" in otu analyses

# TODO compress nosingle.txt and 7_no_primers.fasta.derep

# TODO move all scripts into subdirectories

# TODO add vsearch clustering

# TODO add library-specific size variable

# TODO An attempt to cause the script to exit if any of the commands returns a non-zero status (i.e. FAILS).

# TODO LIBS_ARRAY is never used
# TODO make LIBRARY_DIRECTORIES an array by wrapping it in ()
# TODO if LIBRARY_DIRECTORIES is an array, its length is "${#LIBRARY_DIRECTORIES[@]}"
# TODO for i in "${LIBRARY_DIRECTORIES[@]}"; do echo "${i##*/}" ; done

# TODO LIB_TAG_MOD originally contained sort | uniq; this is unnecessary I think

# TODO add `trap "killall background" EXIT` or `trap 'kill $(jobs -p)' EXIT` to kill background processes on exit

# TODO Add decontamination script

# TODO Add normalization

# TODO streamline config file

# TODO Incorporate warnings for missing columns in metadata

# TODO Add library names variable (require alpha?)

# TODO write table:
 # - rows: libraries, tags (including rows for whole libraries)
 # - column: numbers for number of sequences: successfully merged, filtered, forward index, both indexes (parallelizable using grep during awk demultiplexing?), primer 1, primer 2, singletons, (dups? otus?)

# TODO Find primer index as a function of 6bp preceeding or following primer seq  `grep -E -o '^.{6}'"${primer_F}"''`

# TODO add --debug / --verbose flag to generate files.
# TODO Potential files to be removed as part of cleanup at the end of script:
# - homopolymer_line_numbers.txt
# - 2_filtered.fasta
# - 1_merged.assembled.fastq
# - 1_merged.assembled_A.fastq
# - 1_merged.assembled_B.fastq



# TODO remove whitespace from sequence labels?
# sed 's/ /_/'


# Copy sequences to fasta files into separate directories based on tag sequence on left side of read
# TODO test for speed against removing the tag while finding it: wrap first tag regex in gsub(/pattern/,""):  awk 'gsub(/^.{0,9}'"$TAG_SEQ"'/,""){if . . .


# TODO add single if/else for CONCATENATE_SAMPLES: assign directory as appropriate, correct references within loop to be extendable
if [ "$CONCATENATE_SAMPLES" = "YES" ]; then
	# do the stuff here
	WORKING_DIR="${CONCAT_DIR}"
else
	WORKING_DIR="${LIBRARY_DIRECTORIES}"
fi
################################################################################
# CONCATENATE SAMPLES
################################################################################
# TODO could move this first step up above any loops (no else)
if [ "$CONCATENATE_SAMPLES" = "YES" ]; then

	# TODO MOVE THE VARIABLE ASSIGNMENT TO TOP; MOVE MKDIR TO TOP OF CONCAT IF LOOP
	echo $(date +%H:%M) "Concatenating fasta files..."
	CONCAT_DIR="$ANALYSIS_DIR"/all_lib
	mkdir "${CONCAT_DIR}"
	CONCAT_FILE="${CONCAT_DIR}"/1_demult_concat.fasta

	# TODO could move this into above loop after demultiplexing?
	for CURRENT_LIB in $LIBRARY_DIRECTORIES; do

		LIB_OUTPUT_DIR="${ANALYSIS_DIR}"/${CURRENT_LIB##*/}

		# TODO !!! This will fail if there are underscores in the library names !!!
		# an attempt at making this robust to underscores
		# grep -E -o '_lib_.+?(?=_tag)_tag_.{6}' "${CONCAT_DIR}"/1_demult_concat.fasta | sed 's/_lib_//;s/_tag_/ /' | sort | uniq -c | sort -nr > "${CONCAT_DIR}"/1_demult_concat.fasta.tags

		# TODO wrap in '( ) &' to force into background and allow parallel processing
		# i.e.
		# for primer in "${primers_arr[@]}"; do
		# 	( cutadapt -g ^"${primer}" -e "${PRIMER_MISMATCH_PROPORTION}" -m "${LENGTH_ROI_HALF}" --discard-untrimmed "${CONCAT_FILE}" > "${CONCAT_DIR}"/5_L"${primer}"_removed.fasta ) &
		# done
		# wait


		# TODO rename preliminary to OTU analyses; move analysis script to OTU analysis directory
		OUTPUT_PDF="${ANALYSIS_DIR}"/analysis_results_"${START_TIME}".pdf



# TODO check for all dependencies
if hash pigz 2>/dev/null; then
	ZIPPER="pigz"
	echo "pigz installation found"
else
	ZIPPER="gzip"
	echo "pigz installation not found; using gzip"
fi






exit
# TODO fix this:
# Failed to parse command: export what = DSV format = readname_taxonname separator = comma file = / Users / threeprime / Desktop / Analysis_20151017_1752 / all_lib / meganout_Genus.csv
# Command: quit;
# Executing: quit;
# /Applications/banzai/banzai.sh: line 918: /Users/threeprime/Desktop/Analysis_20151017_1752/all_lib/meganout_Family.csv: No such file or directory
# /Applications/banzai/banzai.sh: line 919: /Users/threeprime/Desktop/Analysis_20151017_1752/all_lib/meganout_Genus.csv: No such file or directory
# Error in file(file, "rt") : cannot open the connection
# Calls: read.csv -> read.table -> file
# In addition: Warning message:
# In file(file, "rt") : cannot open file 'NA': No such file or directory
# Execution halted
# 17:54 passing args to R for preliminary analysis...
# Error in plot.window(xlim, ylim, log, ...) : need finite 'ylim' values
# Calls: stripchart -> stripchart.default -> plot.window
# In addition: Warning messages:
# 1: In min(x) : no non-missing arguments to min; returning Inf
# 2: In max(x) : no non-missing arguments to max; returning -Inf
# Execution halted
# There was a problem generating the PDF.

# TODO documentation:
# shamelessly stolen from https://github.com/geraldinepascal/FROGS
Legend for the next schemas:
    .: Complete nucleic sequence
    !: Region of interest
    *: PCR primers

Paired-end classical protocol:
    In the paired-end protocol R1 and R2 must share a nucleic region.
    For example the amplicons on V3-V4 regions can have a length between
    350 and 500nt; with 2*300pb sequencing the overlap is between 250nt
    and 100nt.
    From:                                    To:
     rDNA .........!!!!!!................    ......!!!!!!!!!!!!!!!!!!!.....
     Ampl      ****!!!!!!****                  ****!!!!!!!!!!!!!!!!!!!****
       R1      --------------                  --------------
       R2      --------------                               --------------

    The maximum overlap between R1 and R2 can be the complete overlap.
        Inconvenient maximum overlap:
        R1    --------------
        R2   --------------
    In this case it is necessary to trim R1 and R2 ends before the process.

    The minimum overlap between R1 and R2 can have 15nt. With less the
    overlap can be incorrect.

Single-end classical protocol:
    rDNA .........!!!!!!................
    Ampl      ****!!!!!!****
    Read      --------------

Custom protocol
    rDNA .....!!!!!!!!!!!!!!............
    Ampl      ****!!!!!!****
    Read      --------------

Note: The amplicons can have a length variability.
      The R1 and R2 can have different length.
