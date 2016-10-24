#!/usr/bin/env bash

# given a large (1.2GB, 8803948 sequences) file:

infile="/Users/threeprime/7_no_primers.fasta"

# one hundred
head -n 100 $infile > temp.fasta
vsearch --derep_fulllength temp.fasta --fasta_width 0 --sizeout --output temp_derep.fasta
# works fine
head temp_derep.fasta

# one thousand
head -n 1000 $infile > temp.fasta
vsearch --derep_fulllength temp.fasta --fasta_width 0 --sizeout --output temp_derep.fasta
head temp_derep.fasta

# one million
head -n 1000000 $infile > temp.fasta
vsearch --derep_fulllength temp.fasta --fasta_width 0 --sizeout --output temp_derep.fasta
head temp_derep.fasta

# seven million
head -n 7000000 $infile > temp.fasta
vsearch --derep_fulllength temp.fasta --fasta_width 0 --sizeout --output temp_derep.fasta
head temp_derep.fasta

# full
echo $infile > temp.fasta
vsearch --derep_fulllength temp.fasta --fasta_width 0 --sizeout --output temp_derep.fasta
head temp_derep.fasta
# now it works?!

