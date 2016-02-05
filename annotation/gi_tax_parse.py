#!/usr/bin/env python

# LOAD SOME MODULES
import gzip

import pandas

# SET SOME FILE PATHS
# for a gi_file
my_gi_file = '/Users/jimmy.odonnell/NCBI/databases/taxonomy/gi_joe.txt'

# or tabular blast output
blast_tab_file = "/Users/jimmy.odonnell/Downloads/blast_20151204/full_output.txt"

# and a gi_taxid file (tab separated, gi number in column 1, taxid in column 2)
gi_taxid_master_file = '/Users/jimmy.odonnell/NCBI/databases/taxonomy/gi_taxid_10000.dmp.gz'


# DO SOME STUFF
# from blast tabular:
my_table = pandas.read_table(blast_tab_file, sep='\t')
gi_column = my_table.ix[:,1]
gi_col_split = gi_column.str.split('|').tolist()
my_gi = set([x[1] for x in gi_col_split])


# from file with 1 gi per line:
my_gi = set(line.strip() for line in open(my_gi_file))

gi_taxid_output_file = '{0:s}.gitax'.format(my_gi_file)

gi_taxid_output = open(gi_taxid_output_file, 'w')


with gzip.open(gi_taxid_master_file, 'r') as gi_taxid:
    for line in gi_taxid:
        line_split = line.rstrip("\n").split('\t')
        if str(line_split[0]) in my_gi:
            print >> gi_taxid_output, (line_split[0] + "\t" + line_split[1])

exit()
# old stuff:
    # found = 0
        # print line
            # print line_split
        # print >> gi_taxid_output, line_split
        # # print line_split[0], found
        #     found = found + 1
        #     print >> gi_taxid_output, line
