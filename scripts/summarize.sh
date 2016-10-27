#!/usr/bin/env bash

banzai_logfile="${1}"

# total reads
total_reads=$( \
  grep -F 'Assembled reads .' "${banzai_logfile}" |\
  awk '{ print $6 }' |\
  sed 's/,//g' |\
  awk '{sum+=$1} END {print sum}'
)

# assembled reads
assembled_reads=$(\
  grep -F 'Assembled reads .' "${banzai_logfile}" |\
  awk '{ print $4 }' |\
  sed 's/,//g' |\
  awk '{sum+=$1} END {print sum}'
)

assembled_perc=$( \
  bc <<< " scale=10; ("${assembled_reads}" / "${total_reads}")*100 "
)

# primary indexing

# expected error filtering
passed_ee_filter=$( \
  fgrep 'sequences kept (of which ' "${banzai_logfile}" |\
  awk '{sum+=$1} END {print sum}'
)
passed_ee_perc=$( \
  bc <<< " scale=10; ("${passed_ee_filter}" / "${assembled_reads}" ) * 100 "
)

# primer index demultiplexing
demult_primers=$( \
  fgrep 'sequences found in primer removal input' "${banzai_logfile}" |\
  awk '{print $3}'
)

demult_perc=$( \
  bc <<< "scale=10; ( "${demult_primers}" / "${passed_ee_filter}" ) * 100 "
)

# primer removal
primers_found=$( \
  awk '/Trimmed reads:/ { print $3 }' "${banzai_logfile}" |\
  tail -n 2 | \
  awk '{sum+=$0} END {print sum}'
)

primers_perc=$( \
  bc <<< "scale=10; ( "${primers_found}" / "${demult_primers}" ) * 100 "
)

# no singletons
no_singletons=$( \
  grep -A 5 'Detecting chimeras' $banzai_logfile | \
  tail -n 1 | \
  awk '{ print $7}'
)

non_chimeras=$( \
  grep -A 5 'Detecting chimeras' $banzai_logfile | \
  tail -n 2 | \
  grep 'non-chimeras' | \
  awk '{ print $(NF-2), $(NF-1)}'
)

# unique sequences
unique_seqs=$( \
  grep -A 2 'Detecting chimeras' $banzai_logfile | \
  grep 'non-chimeras' | \
  awk '{ print $(NF-2)}'
)

echo "A total of" "${total_reads}" "reads in each direction passed the Illumina quality filter and could be assigned to the correct sample on the basis of their primary index."
echo "Of these, ""${assembled_reads}"" ("${assembled_perc:0:5}"%) met the criteria for paired-end assembly and initial quality filtering."
echo "library_reads" "library_percent" "reads were assigned to one of the on the basis of the index sequences ligated during library preparation."
echo "${passed_ee_filter}"" ("${passed_ee_perc:0:5}"%)" "of the assembled reads passed filtering based on expected errors."
echo "${demult_primers}" "("${demult_perc:0:5}"%)" "of the filtered reads contained secondary (primer) index sequences in the expected positions, and "
echo "${primers_found}" "("${primers_perc:0:5}"%)" "of the demultiplexed reads contained primer sequences in the expected position"
echo "Of these reads, a total of" "${no_singletons}" "sequences occurred more than once, "${non_chimeras}" of which passed chimera filtering, in total comprising" "${unique_seqs}" "unique sequences."
