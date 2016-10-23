#!/usr/bin/env bash

# check for chimeras using vsearch

# input: a fasta file. Could be dereplicated or OTU clustered sequences. Most folks seem to think its best to check for chimeras after OTU clustering.
# input fasta file should present abundance annotations in the fasta header pattern ([;]size=integer[;])
# usage: bash chimera_check.sh OTUs.fasta

infile="${1}"
chimera_free_fasta="${infile%.*}"_nochime.fasta
chimera_fasta="${infile%.*}"_chime.fasta
borderline_fasta="${infile%.*}"_borderchime.fasta
uchime_outfile="${infile%.*}".uchime
align_outfile="${infile%.*}".aln

vsearch \
  --uchime_denovo "${infile}" \
  --nonchimeras "${chimera_free_fasta}" \
  --chimeras "${chimera_fasta}" \
  --uchimealns "${align_outfile}" \
  --uchimeout "${uchime_outfile}" \
  --borderline "${borderline_fasta}" \
  --fasta_width 0 \
  --alignwidth 0 \
  --sizeout \
  --abskew 2.0 \
  --mindiffs 3 \
  --mindiv 0.8 \
  --minh 0.28


# --uchimeout filename
# Write chimera detection results to filename using the uchime tab-separated format of 18 fields (see the list below). Use --uchimeout5 to use a format compatible with use- arch v5 and earlier versions. Rows output order may vary when using multiple threads.
# 1. score: higher score means a more likely chimeric alignment.
# 2. Q: query sequence label.
# 3. A: parent A sequence label.
# 4. B: parent B sequence label.
# 5. T: top parent sequence label (i.e. parent most similar to the query). That field is removed when using --uchimeout5.
# 6. idQM: percentage of similarity of query (Q) and model (M) constructed as a part of parent A and a part of parent B.
# 7. idQA: percentage of similarity of query (Q) and parent A.
# 8. idQB: percentage of similarity of query (Q) and parent B.
# 9. idAB: percentage of similarity of parent A and parent B.
# 10. idQT: percentage of similarity of query (Q) and top parent (T).
# 11. LY: yes votes in the left part of the model.
# 12. LN: no votes in the left part of the model.
# 13. LA: abstain votes in the left part of the model.
# 14. RY: yes votes in the right part of the model.
# 15. RN: no votes in the right part of the model.
# 16. RA: abstain votes in the right part of the model.
# 17. div: divergence, defined as (idQM - idQT).
# 18. YN: query is chimeric (Y), or not (N), or is a borderline case (?).
