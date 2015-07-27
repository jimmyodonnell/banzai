#banzai!#

**banzai** is a BASH shell script that links together the disparate analyses needed to process the raw sequencing results from an Illumina run, the eventual goal being a table of the number of sequences per taxon found in a set of samples. Preliminary ecological analyses are included as well.

The script should run on Unix and Linux machines. The script makes heavy usage of Unix command line utilities (such as find, grep, sed, awk, and more) and is written for the BSD versions of those programs as found on standard installations of Mac OSX. I tried to use POSIX-compliant commands wherever possible.

## Basic implementation ##
Simply edit the parameters in the file 'banzai_params.sh', then type into a terminal:

```sh
bash /Users/user_name/path/to/the/file/banzai.sh
```

It's important to use `bash` rather than `sh` or `.` to invoke the script. Someday I'll figure out a better workaround, but for now this was the only way I could guarantee the log file was created in the way I wanted.

### This could take a while... ###
In Mac OS Mountain Lion and later, you can override your computer's sleep settings by running the script like so:

```sh
caffeinate -i -s bash /Users/user_name/path/to/the/file/banzai.sh
```

## Dependencies ##
Aside from the standard command line utilities that are already included on Unix machines (awk, sed, grep, etc), this script relies on the following tools:

* **PEAR**: merging paired-end reads
* **usearch**: quality filtering of merged paired-end fastq files, OTU clustering (could be replaced with vsearch and/or swarm)
* **cutadapt**: primer removal (I might replace with awk)
* **seqtk**: reverse complementing entire fastq/a files
* **python**: fast consolidation of duplicate sequences (installed by default on Macs)
* **blast+**: taxonomic assignment
* **MEGAN**: taxonomic assignment
* **R**: ecological analyses. Requires the packages **vegan** and **gtools**

I might also consider adding:

* **fastqc**: quality control of raw sequencing fastq files


## Sequencing Pool Metadata ##
If you provide a CSV spreadsheet that contains metadata about the samples, banzai can read some of the parameters from it, like the primers and multiplex index sequences. You need to provide the file path to the spreadsheet, and the relevant column names.

It is VERY important that this file be encoded with UNIX line breaks. You can do this from Excel and TextWrangler. It doesn't appear to be critical that the text is encoded using UTF-8, though this is certainly the safest option. Early in the logfile you can check to be sure the correct number of tags and primer sequences were found.


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

Note that the original script also ouput a file of the sequences only (no names), but I removed this functionality on 20150417

##Wishlist/TODO/notes to self##
* Add decontamination script
* Add normalization
* streamline config file
* Incorporate warnings for missing columns in metadata
* Add library names variable (require alpha?)
* write table:
 - rows: libraries, tags (including rows for whole libraries)
 - column: numbers for number of sequences: successfully merged, filtered, forward index, both indexes (parallelizable using grep during awk demultiplexing?), primer 1, primer 2, singletons, (dups? otus?)

* Find primer index as a function of 6bp preceeding or following primer seq  `grep -E -o '^.{6}'"${primer_F}"''`

* Potential files to be removed as part of cleanup at the end of script:
 - homopolymer_line_numbers.txt
 - 2_filtered.fasta
 - 1_merged.assembled.fastq
 - 1_merged.assembled_A.fastq
 - 1_merged.assembled_B.fastq

## Known Issues/Bugs ##
* As of 20150614, libraries must be in folders called lib1, lib2, etc. Need to fix the sed renaming scheme to accommodate variable library names

###Notes###

2014 11 12: Noticed that the reverse tag removal step removes the tag label from the sequenceID line of fasta files if the tag sequence is RC-palindromic!

An alternate hack to have the pipeline print to terminal AND file, in case logging breaks:
sh script.sh  2>&1 | tee ~/Desktop/logfile.txt
