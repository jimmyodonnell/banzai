#banzai!#

🏄

**banzai** is a BASH (shell) script that links together the disparate programs needed to process the raw sequencing results from an Illumina run into a contingency table of the number of sequences per taxon found in a set of samples. Some preliminary ecological analyses are included as well.

The script should run on Unix and Linux machines. The script makes heavy usage of Unix command line utilities (such as find, grep, sed, awk, and more) and is written for the BSD versions of those programs as found on standard installations of Mac OSX. I tried to use POSIX-compliant commands wherever possible.

## Basic implementation ##
**NEW!!!**
*NOTE* that as of 2015-10-09, you must direct banzai.sh to your parameter file. This allows for much easier use when analyzing multiple types of projects. Parameter files can be called whatever you want -- e.g. `banzai_params_16s.sh`. When you invoke the file banzai.sh, it will source whatever file you give it using the first argument (separated by a space).
Simply copy the file 'banzai_params.sh' into a new folder, set parameters as desired, then type into a terminal:

```sh
bash /Users/user_name/path/to/the/file/banzai.sh   /User/user_name/path/to/param_file.sh
```

It's important to use `bash` rather than `sh` or `.` to invoke the script. Someday I'll figure out a better workaround, but for now this was the only way I could guarantee the log file was created in the way I wanted.


## Dependencies ##
Aside from the standard command line utilities (awk, sed, grep, etc) that are already included on Unix machines, this script relies on the following tools:

* **[PEAR](http://sco.h-its.org/exelixis/web/software/pear/)**: merging paired-end reads
* **[cutadapt](https://github.com/marcelm/cutadapt)**: primer removal (I might replace with awk)
* **[vsearch](https://github.com/torognes/vsearch)**: sequence quality filtering (requires version 1.4.0 or greater); OTU clustering
* **[swarm](https://github.com/torognes/swarm)**: OTU clustering
* **[seqtk](https://github.com/lh3/seqtk)**: reverse complementing entire fastq/a files
* **[python](https://www.python.org/)**: fast consolidation of duplicate sequences (installed by default on Macs)
* **[blast+](http://www.ncbi.nlm.nih.gov/books/NBK279690/)**: taxonomic assignment
* **[MEGAN](http://ab.inf.uni-tuebingen.de/software/megan5/)**: taxonomic assignment
* **[R](https://www.r-project.org/)**: ecological analyses. Requires the packages **vegan** and **gtools**

### Recommended ###
* Compressing and decompressing files can be slow because standard, built-in utilities (gzip) do not run in parallel. Installing the parallel compression tool **[pigz](http://zlib.net/pigz/)** can yield substantial speedups. Banzai will check for pigz and use it if available.

* I recommend that before analyzing data, you check and report basic properties of the sequencing runs using **[fastqc](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)**. I have included a script to do this for all the fastq or fastq.gz files in any subdirectory of a directory (run_fastqc.sh).

### Optional/Deprecated ###
* **[usearch](http://www.drive5.com/usearch/)**: filtering paired reads on the basis of the sum of the error probabilities (maximum expected errors). This can be turned off, probably without much change in final data quality. We used to do OTU clustering with usearch, but the 32bit version can't handle larger data sets.

## Sequencing Pool Metadata ##
If you provide a CSV spreadsheet that contains metadata about the samples, banzai can read some of the parameters from it, like the primers and multiplex index sequences. You need to provide the file path to the spreadsheet, and the relevant column names.

It is VERY important that this file be encoded with UNIX line breaks. You can do this from Excel and TextWrangler. It doesn't appear to be critical that the text is encoded using UTF-8, though this is certainly the safest option. Early in the logfile you can check to be sure the correct number of tags and primer sequences were found.

No field should contain any spaces. That means row names, column names, and cells. Accomodating this would require an advanced degree in bash-quoting judo, which I do not have.

## LIBRARY NAMES ##
As of 2015-10-09, libraries no longer have to be named anything in particular (e.g. A, B, lib1, lib2),
BUT THEY CANNOT CONTAIN UNDERSCORES or spaces!

## Organization of raw data ##
Your data (fastq files) can be compressed or not; but banzai currently only works with paired-end Illumina data. Thus, the bare minimum input is two fastq files corresponding to the first and second read. *Banzai will fail if there are files in your library folders that are not your raw data but have 'fastq' in the filename!* For example, if your library contains four files: "R1.fastq", "R1.fastq.gz", "R2.fastq", and "R2.fastq.gz". banzai will grab the first two (R1.fastq and R1.fastq.gz) and try to merge them, and (correctly) fail miserably. Note that while PEAR 0.9.7 merges compressed (\*.gz) files directly, PEAR 0.9.6 does not do so correctly. If given compressed files as input, banzai first decompresses them, which will add a little bit of time to the overall analysis.

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
* Find primer index as a function of 6bp preceeding or following primer seq  `grep -E -o '^.{6}'"${primer_F}"''`
* Potential files to be removed as part of cleanup at the end of script:
 - homopolymer_line_numbers.txt
 - 2_filtered.fasta
 - 1_merged.assembled.fastq
 - 1_merged.assembled_A.fastq
 - 1_merged.assembled_B.fastq



### This could take a while... ###
In Mac OS Mountain Lion and later, you can override your computer's sleep settings by running the script like so:

```sh
caffeinate -i -s bash /Users/user_name/path/to/the/file/banzai.sh
```

## Known Issues/Bugs ##
* As of 20150614, libraries must be in folders called lib1, lib2, etc. Need to fix the sed renaming scheme to accommodate variable library names

###Notes###

2014 11 12: Noticed that the reverse tag removal step removes the tag label from the sequenceID line of fasta files if the tag sequence is RC-palindromic!

An alternate hack to have the pipeline print to terminal AND file, in case logging breaks:
sh script.sh  2>&1 | tee ~/Desktop/logfile.txt

* 2015-10-09 read length calculated from raw data. Library names are flexible.
