# banzai!

🏄

**banzai** is a shell script (bash) that links together the disparate programs needed to process the raw results from an Illumina sequencing run of PCR amplicons into a contingency table of the number of similar sequences found in each of a set of samples.

## Contents
1. [Introduction](#introduction)
2. [Flow Chart](#flow-chart)
3. [Usage](#usage)
4. [Dependencies](#dependencies)
5. [Sequencing Metadata](#sequencing-metadata)
5. [Bugs and Issues](#bugs-and-issues)
6. [Notes](#notes)

## Introduction
The script should run on Unix (Mac OSX) and Linux machines. It makes heavy usage of Unix command line utilities (such as find, grep, sed, awk, and more) and is written for the BSD versions of those programs as found on standard installations of Mac OSX. I tried to use POSIX-compliant commands wherever possible.

### Lab Preparation: Samples to sequences ###
Banzai was designed for sequencing data that were generated for a target region like so:

```
Genomic DNA:    ----------------------------------------------------------------
target region:                   ~~~~~~~~~~~~~~~~~~~~~~~~
```

**PCR:**
```
Primers:                   ******                        ******
Secondary index:        +++                                    +++
full primer:            +++******                        ******+++
amplicon:               +++******~~~~~~~~~~~~~~~~~~~~~~~~******+++
```
The amplicon contains everything that will show up in the sequence. This is more generally referred to as the 'insert'. Note that the insert size is greater than the size of the region of interest alone.

**Library Prep:**
```
primary index:       :::                                          :::
adapter:           aa                                                aa
final fragment:    aa:::+++******~~~~~~~~~~~~~~~~~~~~~~~~******+++:::aa
```
This is the fragment that actually goes into the sequencer. Note that it is bigger than the insert (above), and that the adapter and primary indexes shouldn't end up in your sequences.

**Sequencing:**
```
Read 1:            aa:::+++******~~~~~~~~~~~~~~
Read 2:                                    ~~~~~~~~~~~~~~******+++:::aa
```

### Banzai: Sequences to Samples ###
Banzai works backward through the same process:

**Demultiplexing (primary):**
```
Read 1:                 +++******~~~~~~~~~~~~~~
Read 2:                                    ~~~~~~~~~~~~~~******+++
```
(This has probably already been performed by the time you get the data.)

**Read merging:**
```
merged reads:           +++******~~~~~~~~~~~~~~~~~~~~~~~~******+++
```

**Demultiplexing (secondary):**
```
                           ******~~~~~~~~~~~~~~~~~~~~~~~~******
 ```

**Primer removal:**
```
                                 ~~~~~~~~~~~~~~~~~~~~~~~~
```
*(Layout inspired in part by the [FROGS](https://github.com/geraldinepascal/FROGS) documentation.)*


## Flow Chart ##
![](doc/banzai_flowchart.png)

## Usage ##
Copy the file 'banzai_params.sh' into a new folder and set parameters as desired. Then run the banzai script, using your newly edited parameter file like so (Mac OSX):

```sh
bash /Users/user_name/path/to/the/file/banzai.sh   /User/user_name/path/to/param_file.sh
```

It's probably important to use `bash` rather than `sh` or `.` to invoke the script. Someday I'll figure out a better workaround, but for now this was the only way I could guarantee the log file was created in the way I wanted.


## Dependencies ##
Aside from the standard command line utilities (awk, sed, grep, etc) that are already included on Unix machines, this script relies on the following tools:

* **[PEAR](http://sco.h-its.org/exelixis/web/software/pear/)**: merging paired-end reads
* **[cutadapt](https://github.com/marcelm/cutadapt)**: primer removal (I might replace with awk)
* **[vsearch](https://github.com/torognes/vsearch)**: sequence quality filtering (requires version 1.4.0 or greater); OTU clustering
* **[swarm](https://github.com/torognes/swarm)**: OTU clustering
* **[seqtk](https://github.com/lh3/seqtk)**: reverse complementing entire fastq/a files
* **[python](https://www.python.org/)**: fast consolidation of duplicate sequences (installed by default on Macs)
* **[blast+](http://www.ncbi.nlm.nih.gov/books/NBK279690/)**: taxonomic assignment
* **[R](https://www.r-project.org/)**: ecological analyses. Requires these packages:
  * [data.table](https://cran.r-project.org/web/packages/data.table/index.html): data manipulation. Default installation on Mac probably won't enable parallel capabilities; to do so, check [these instructions](https://github.com/Rdatatable/data.table/wiki/Installation).
  * [gtools](https://cran.r-project.org/web/packages/gtools/index.html): data manipulation
  * [reshape2](https://cran.r-project.org/web/packages/reshape2/index.html): data manipulation
  * [vegan](https://cran.r-project.org/web/packages/vegan/index.html): ecological analyses
  * [taxize](https://cran.r-project.org/web/packages/taxize/index.html): taxonomic annotation


Follow the [Vagrant-VirtualBox instructions](doc/vagrant_install.md) to automatically install your own virtual machine that includes all of these dependencies.

### Recommended ###
* Compressing and decompressing files can be slow because standard, built-in utilities (gzip) do not run in parallel. Installing the parallel compression tool **[pigz](http://zlib.net/pigz/)** can yield substantial speedups. Banzai will check for pigz and use it if available.

* I recommend that before analyzing data, you check and report basic properties of the sequencing runs using **[fastqc](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)**. I have included a script to do this for all the fastq or fastq.gz files in any subdirectory of a directory (run_fastqc.sh).


## Sequencing Metadata ##
This is a critical component. 
You must provide a CSV spreadsheet that contains metadata about the samples. Banzai will read some of the parameters from it, like the primers and multiplex index sequences. You need to provide the file path to the spreadsheet, and the relevant column names. Some of the details seem tedious (like listing each file name), but they are inspired by the EMBL/EBIMetagenomics metadata requirements. That is, you're going to have to do this stuff at some point anyway.

This file should be encoded with UNIX newline characters (LF). Banzai will attempt to check for and fix files encoded with Windows newline characters (CRLF), but this is a place to look if you get mysterious errors. Early in the logfile you can check to be sure the correct number of tags and primer sequences were found.

No field should contain any spaces. That means row names, column names, and cells. Accommodating this would require an advanced degree in bash-quoting judo, which I do not have.


## Bugs and Issues ##
* Currently awaiting catastrophic finding...
- be aware of potential problems between awk versions on Linux and Mac. On a Mac, the output of `awk --version` is `awk version 20070501`


## Notes ##
An alternate hack to have the pipeline print to terminal AND file, in case logging breaks:
sh script.sh  2>&1 | tee ~/Desktop/logfile.txt

* 2016-10-22 Began major reorganization. Created v0.1.0-beta for MBON in case anything causes a fire.
* 2015-10-19 expected error filtering implemented via vsearch. OTU clustering can be done with swarm or usearch.
* 2015-10-09 read length calculated from raw data. Library names are flexible.
* 2014-11-12 Noticed that the reverse tag removal step removed the tag label from the sequenceID line of fasta files if the tag sequence is RC-palindromic!
* **Library Names**:
As of 2015-10-09, libraries no longer have to be named anything in particular (e.g. A, B, lib1, lib2),
BUT THEY CANNOT CONTAIN UNDERSCORES or spaces! (This will be moot once library index sequences are required)
* **Raw Data**: 
Your data (fastq files) can be compressed or not; but banzai currently only works with paired-end Illumina data. 
Thus, the bare minimum input is two fastq files corresponding to the first and second read. *Banzai will fail if there are files in your library folders that are not your raw data but have 'fastq' in the filename!* 
For example, if your library contains four files: "R1.fastq", "R1.fastq.gz", "R2.fastq", and "R2.fastq.gz". banzai will grab the first two (R1.fastq and R1.fastq.gz) and try to merge them, and (correctly) fail miserably. Note that while PEAR 0.9.7 merges compressed (\*.gz) files directly, PEAR 0.9.6 does not do so correctly. If given compressed files as input, banzai first decompresses them, which will add a little bit of time to the overall analysis.
