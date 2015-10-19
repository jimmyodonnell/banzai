#!/usr/bin/env bash



# TODO add case statement for clustering
cluster_method="swarm" #[ swarm | vsearch | usearch ]

case "${cluster_method}" in

    "swarm" )

        source "${SCRIPT_DIR}"/OTU_clustering/cluster_swarm.sh

    ;;

    "vsearch" )

        source "${SCRIPT_DIR}"/OTU_clustering/cluster_vsearch.sh

    ;;

    "usearch" )

        source "${SCRIPT_DIR}"/OTU_clustering/cluster_usearch.sh

    ;;

    * )

        echo "${cluster_method}" 'is an invalid clustering method.'
        echo 'Must be one of swarm, vsearch, or usearch. Using swarm.'
        source "${SCRIPT_DIR}"/OTU_clustering/cluster_swarm.sh

    ;;

esac


# after usearch, check if there was an error:

usearch_exit_status=$( echo $? )
if [[ "${usearch_exit_status}" = 1 ]]; then
  echo 'usearch exited with an error -- something went wrong!'
  echo 'Aborting script.'
  exit
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
