#!/usr/bin/env bash

# check the status of your blast database

EXPORT="FALSE"
if [[ "${EXPORT}" = "TRUE" ]]; then
  # to write to file:
  CURRENT_TIME=$(date +%Y%m%d_%H%M)
  BLAST_REPORT_OUT="blastdb_report_""${CURRENT_TIME}".txt
  blastdbcmd -db nt -info  > "${BLAST_REPORT_OUT}"
else
  blastdbcmd -db nt -info
fi

UPDATE="FALSE"
if [[ "${UPDATE}" = "TRUE" ]]; then

  # find the update script
  DIR=$( dirname $( which blastn ))

  cd $BLASTDB
  perl "${DIR}"/update_blastdb.pl --decompress nt

fi
