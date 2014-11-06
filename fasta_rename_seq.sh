#!/bin/bash

# Rename Illumina-generated fasta sequences.

INFILE='filt_80K.fasta'

sed -E "s/ (1|2):N:0:1/_/" "${INFILE}" > tmp.fasta

sed -E "s/>([a-zA-Z0-9-]*:){4}/>/" tmp.fasta > renamed.fasta

rm tmp.fasta
