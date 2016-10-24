## A note on removal of duplicate sequences##

###  (dereplicate_fasta.py) ###

* Input: a fasta file (e.g. 'infile.fasta')

* Output: a file with the same name as the input but with the added extension '.derep' (e.g. 'infile.fasta.derep')

This output file contains each unique DNA sequence from the fasta file, followed by the labels of the reads matching this sequence
Thus, if an input fasta file consisted of three reads with identical DNA sequences:

	>READ1
	AATAGCGCTACGT
	>READ2
	AATAGCGCTACGT
	>READ3
	AATAGCGCTACGT

The output file is as follows:

	AATAGCGCTACGT; READ1; READ2; READ3

Note that the original script also output a file of the sequences only (no names), but I removed this functionality on 20150417
