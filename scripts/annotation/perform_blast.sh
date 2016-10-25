#!/usr/bin/env bash

################################################################################
# TAXONOMIC ANNOTATION
################################################################################
if [ "$PERFORM_BLAST" = "YES" ]; then
	echo $(date +%Y-%m-%d\ %H:%M) "BLASTing..."
	blast_output="${DEREP_INPUT%/*}"/10_BLASTed.xml
	blastn \
		-query "${BLAST_INPUT}" \
		-db "$BLAST_DB" \
		-num_threads "$n_cores" \
		-perc_identity "${PERCENT_IDENTITY}" \
		-word_size "${WORD_SIZE}" \
		-evalue "${EVALUE}" \
		-max_target_seqs "${MAXIMUM_MATCHES}" \
		-culling_limit "${culling_limit}" \
		-outfmt 5 \
		-out "${blast_output}"

	# check for blast output
	if [[ ! -s "${blast_output}"  ]]; then
		echo
		echo 'BLAST failed: the output file is empty or absent.'
	    echo 'File should be:' "${blast_output}"
		echo
	fi
else
	echo 'BLAST search not performed.'
	echo
fi
