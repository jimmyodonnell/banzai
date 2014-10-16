#!/bin/bash

#my_dir="$(dirname "$0")"

#Standard (mandatory) options:

#INFILE1=/Users/threeprime/Documents/Projects/Puget_Sound_eDNA/Analysis/eDNA_analysis/co1_data_20140716/Raw_Sequence_Data/MiSeq-Stanford_SFGF_July_2014_CO1-Leray_Test/CO1-Tilapia-Tissue/CO1-Tilap_S8_L001_R1_001.fastq
#INFILE1=$(find $my_dir -name '*R1*fastq')
#  -f, --forward-fastq         <str>     Forward paired-end FASTQ file.

#INFILE2=/Users/threeprime/Documents/Projects/Puget_Sound_eDNA/Analysis/eDNA_analysis/co1_data_20140716/Raw_Sequence_Data/MiSeq-Stanford_SFGF_July_2014_CO1-Leray_Test/CO1-Tilapia-Tissue/CO1-Tilap_S8_L001_R2_001.fastq
#INFILE2=$(find $my_dir -name '*R2*fastq')
#  -r, --reverse-fastq         <str>     Reverse paired-end FASTQ file.

#  -o, --output                <str>     Output filename.
# OUTFILE=$my_dir/merged

#Optional:

OVERLAP_EXPECTED=$(($LENGTH_FRAG - (2 * ($LENGTH_FRAG - $LENGTH_READ) ) ))
# Alternate: OVERLAP_EXPECTED=$(($LENGTH_READ - ($LENGTH_FRAG - $LENGTH_READ)))

MINOVERLAP=$(( $OVERLAP_EXPECTED / 2 ))
#  -v, --min-overlap           <int>     Specify the minimum overlap size. The minimum overlap may be
#                                        set to 1 when the statistical test is used. However, further
#                                        restricting  the  minimum overlap size to a proper value may
#                                        reduce false-positive assembles. (default: 10)

ASSMAX=$(( $LENGTH_FRAG + 10 ))
#  -m, --max-assembly-length   <int>     Specify   the  maximum  possible  length  of  the  assembled
#                                        sequences.  Setting this value to 0 disables the restriction
#                                        and assembled sequences may be arbitrary long. (default: 0)

ASSMIN=$(( $LENGTH_FRAG - 10 ))
#  -n, --min-assembly-length   <int>     Specify   the  minimum  possible  length  of  the  assembled
#                                        sequences.  Setting this value to 0 disables the restriction
#                                        and  assembled  sequences  may be arbitrary short. (default:
#                                        50)

QT=15
#  -q, --quality-threshold     <int>     Specify  the  quality  score  threshold for trimming the low
#                                        quality  part  of  a  read.  If  the  quality  scores of two
#                                        consecutive  bases  are  strictly  less  than  the specified
#                                        threshold,  the  rest of the read will be trimmed. (default:
#                                        0)

TRIMMIN=1
#  -t, --min-trim-length       <int>     Specify  the  minimum length of reads after trimming the low
#                                        quality part (see option -q). (default: 1)

UNCALLEDMAX=1
#  -u, --max-uncalled-base     <float>   Specify  the maximal proportion of uncalled bases in a read.
#                                        Setting this value to 0 will cause PEAR to discard all reads
#                                        containing  uncalled  bases.  The other extreme setting is 1
#                                        which  causes  PEAR  to process all reads independent on the
#                                        number of uncalled bases. (default: 1)

TEST=1
#  -g, --test-method           <int>     Specify  the  type  of  statistical  test.  Two  options are
#                                        available. (default: 1)
#                                        1: Given the minimum allowed overlap, test using the highest
#                                        OES. Note that due to its discrete nature, this test usually
#                                        yields  a lower p-value for the assembled read than the cut-
#                                        off  (specified  by -p). For example, setting the cut-off to
#                                        0.05  using  this  test,  the  assembled reads might have an
#                                        actual p-value of 0.02.
#
#                                        2. Use the acceptance probability (m.a.p). This test methods
#                                        computes  the same probability as test method 1. However, it
#                                        assumes  that  the  minimal  overlap is the observed overlap
#                                        with  the  highest  OES, instead of the one specified by -v.
#                                        Therefore,  this  is  not  a  valid statistical test and the
#                                        'p-value'  is  in fact the maximal probability for accepting
#                                        the assembly. Nevertheless, we observed in practice that for
#                                        the case the actual overlap sizes are relatively small, test
#                                        2  can  correctly  assemble  more  reads  with only slightly
#                                        higher false-positive rate.


PVALUE=0.01
#  -p, --p-value               <float>   Specify  a p-value for the statistical test. If the computed
#                                        p-value of a possible assembly exceeds the specified p-value
#                                        then  paired-end  read  will not be assembled. Valid options
#                                        are: 0.0001, 0.001, 0.01, 0.05 and 1.0. Setting 1.0 disables
#                                        the test. (default: 0.01)


#  -e, --empirical-freqs                 Disable  empirical base frequencies. (default: use empirical
#                                        base frequencies)

SCORING=2
#  -s, --score-method          <int>     Specify the scoring method. (default: 2)
#                                        1. OES with +1 for match and -1 for mismatch.
#                                        2: Assembly score (AS). Use +1 for match and -1 for mismatch
#                                        multiplied by base quality scores.
#                                        3: Ignore quality scores and use +1 for a match and -1 for a
#                                        mismatch.

#  -b, --phred-base            <int>     Base PHRED quality score. (default: 33)

# NOTE THAT MEMORY >200MB CAN _______SLOW DOWN__________ THE ANALYSIS!!!
#  -y, --memory                <str>     Specify  the  amount of memory to be used. The number may be
#                                        followed  by  one  of  the  letters  K,  M,  or  G  denoting
#                                        Kilobytes,  Megabytes and Gigabytes, respectively. Bytes are
#                                        assumed in case no letter is specified.

#  -c, --cap                   <int>     Specify  the upper bound for the resulting quality score. If
#                                        set to zero, capping is disabled. (default: 40)

THREADS=16
#  -j, --threads               <int>     Number of threads to use
# NOTE THAT THREADS SHOULD BE SET FOR NUMBER OF CORES ON MACHINE (MULTIPLIED BY HYPERTHREAD CAPACITY).
# CAN BE SET HIGHER WITHOUT SUBSTANTIAL SPEED LOSS.
#  -h, --help                            This help screen.
