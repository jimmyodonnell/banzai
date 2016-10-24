# Installation Notes #

Installing lots of command line tools requires navigating a certain level of command-line ~~bullshittery~~ wizardry ([source](http://www.pgbovine.net/command-line-bullshittery.htm)). Here are some tips intended to help you get up and running on a Mac.


First things first: Mac users should install command-line tools using the following command:  
`xcode-select --install`


[Homebrew](http://brew.sh/)  
Homebrew isn't used by banzai, but installing it will make your life MUCH easier. Homebrew allows you to download and install a large number of tools right from the command line, avoiding many of the hassles. Install Homebrew in one easy go:
`ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"`


[PEAR](http://sco.h-its.org/exelixis/web/software/pear/)  
`brew install homebrew/science/pear`

[vsearch](https://github.com/torognes/vsearch)  
either:
`brew install https://raw.githubusercontent.com/tseemann/homebrew-science/f55e38c5c941f8cf6db68cd576baafc6a1d89b1e/vsearch.rb`
or
Download binary for v1.4.0
`cp /Users/jimmy.odonnell/Applications/install/vsearch-1.4.0-osx-x86_64/bin/vsearch /usr/bin`
`mv /Users/jimmy.odonnell/Applications/install/vsearch-1.4.0-osx-x86_64/man/vsearch.1.gz /usr/share/man/man1/`

[cutadapt](https://github.com/marcelm/cutadapt)  
Mysteriously disappeared from Jimmy's work computer (! AHA ! I think because I used pip, which was installed by Anaconda, which I then moved!), and Homebrew installation (`brew install cutadapt`) appeared broken as of 2015-10-19. To fix it, used `easy_install`, a python package managing tool that may or may not be automatically included in Mac python:  
`easy_install cutadapt`  
`pip install cutadapt`

[seqtk](https://github.com/lh3/seqtk)  
`brew install seqtk`  
or, download the compressed file, unzip, `cd` into it, then  
`make`  
`cp /Users/jimmy.odonnell/Downloads/seqtk-master/seqtk /usr/bin`

[swarm](https://github.com/torognes/swarm#install)  
`cd /Users/jimmy.odonnell/Applications/install/swarm-master/src`  
`make` (requires gcc)  
`cp ../bin/swarm /usr/bin/swarm`  
`cd ../man`  
`gzip -c swarm.1 > swarm.1.gz`  
`mv swarm.1.gz /usr/share/man/man1/`

[blast+](http://www.ncbi.nlm.nih.gov/books/NBK279690/)  
`brew install homebrew/science/blast`

[pigz](http://zlib.net/pigz/)  
Download .tar.gz file. Double-click to decompress it. `cd` it the new folder. Type:  
`make`  
then  
`cp /Users/jimmy.odonnell/Applications/pigz-2.3.3/pigz /usr/bin`


[fastqc](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)  
unpack dmg file into your Applications folder.
create symbolic link:
`ln -s /Applications/FastQC.app/Contents/Resources/Java/fastqc /usr/bin/fastqc`

[MEGAN](http://ab.inf.uni-tuebingen.de/software/megan5/)  
You're on your own here.

[usearch](http://www.drive5.com/usearch/)  
(no longer required)
`cp /Users/jimmy.odonnell/Downloads/usearch8.1.1756_i86osx32 /usr/bin/usearch`

[R](https://www.r-project.org/)  
Self-explanatory.
