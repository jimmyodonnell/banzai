#!/usr/bin/python

import sys
import fasta
import fastq
import primer

MAX_PRIMER_MISMATCHES = 2

FileName = sys.argv[1]
Primer = sys.argv[2]

SeqCount = 0
OutCount = 0
PrimerMismatchCount = 0

PL = len(Primer)

def MatchesPrimer(Seq, Primer):
	return primer.MatchPrefix(Seq, Primer)

def OnRec(Label, Seq, Qual):
	global PL, SeqCount, OutCount, PrimerMismatchCount #LabelPrefix, Barcode, BarcodeMismatchCount, 
	SeqCount += 1
	Diffs = MatchesPrimer(Seq, Primer)
	if Diffs > MAX_PRIMER_MISMATCHES:
		PrimerMismatchCount += 1
		return
	OutCount += 1
	fastq.WriteRec(sys.stdout, Label, Seq[PL:], Qual[PL:]) #NewLabel, 
fastq.ReadRecs(FileName, OnRec)
print >> sys.stderr, "%10u seqs" % SeqCount
print >> sys.stderr, "%10u matched" % OutCount
print >> sys.stderr, "%10u primer mismatches" % PrimerMismatchCount
