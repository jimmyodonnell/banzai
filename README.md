#banzai!#

**banzai** is a BASH shell script that links together the disparate programs needed to process the raw sequencing results from an Illumina run into a table of the number of sequences per taxon found in a set of samples. Some preliminary ecological analyses are included as well.

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
Aside from the standard command line utilities (awk, sed, grep, etc) that are already included on Unix machines, this script relies on the following tools:

* **[PEAR](http://sco.h-its.org/exelixis/web/software/pear/)**: merging paired-end reads
* **[usearch](http://www.drive5.com/usearch/)**: filtering paired reads on the basis of the sum of the error probabilities (maximum expected errors). This can be turned off, probably without much change in final data quality. We used to do OTU clustering with usearch, but the 32bit version can't handle larger data sets.
* **[swarm](https://github.com/torognes/swarm)**: OTU clustering
* **[vsearch](https://github.com/torognes/vsearch)**: OTU clustering
* **[cutadapt](https://github.com/marcelm/cutadapt)**: primer removal (I might replace with awk)
* **[seqtk](https://github.com/lh3/seqtk)**: reverse complementing entire fastq/a files
* **[python](https://www.python.org/)**: fast consolidation of duplicate sequences (installed by default on Macs)
* **[blast+](http://www.ncbi.nlm.nih.gov/books/NBK279690/)**: taxonomic assignment
* **[MEGAN](http://ab.inf.uni-tuebingen.de/software/megan5/)**: taxonomic assignment
* **[R](https://www.r-project.org/)**: ecological analyses. Requires the packages **vegan** and **gtools**

I recommend that before analyzing data, you check and report basic properties of the sequencing runs using **[fastqc](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)**. I have included a script to do this for all the fastq or fastq.gz files in any subdirectory of a directory (run_fastqc.sh).


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
* add `trap "killall background" EXIT` or `trap 'kill $(jobs -p)' EXIT` to kill background processes on exit
* Add decontamination script
* Add normalization
* streamline config file
* Incorporate warnings for missing columns in metadata
* Add library names variable (require alpha?)
* write table:
 - rows: libraries, tags (including rows for whole libraries)
 - column: numbers for number of sequences: successfully merged, filtered, forward index, both indexes (parallelizable using grep during awk demultiplexing?), primer 1, primer 2, singletons, (dups? otus?)
* 
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
