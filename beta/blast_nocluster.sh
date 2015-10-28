#!/bin/bash

cd '/Users/threeprime/Desktop/MiSeq-Stanford_SFGF-Sept_2014_16S-Eelgrass_Puget_Sound/16S-Eelgrass'

blastn -query 3_derep.fasta -db "/Users/threeprime/Desktop/16Smetazoa/16Smetazoa.fasta" -perc_identity 95 -word_size 50 -max_target_seqs 10 -outfmt 5 -out BLAST_noclust.xml
