#!/usr/bin/env python

'''
or keep appending index numbers to files by match_name, then bring them all together at end

slowdown is mostly IO and searching
'''
import sys

try:
	infile = sys.argv[1]
except:
	raise RuntimeError, '\n\n\n\tusage: ./dereplicate_fasta.py <filename_to_sort>\n\n'

kwargs = {
		#_file to sort
		'fname'		: infile,

		#_what you want between string1 and uniq_id1, uniq_id2 sections
		'sep'		: ' ; ',

		#_name of output file for col1 == dna str, col2-coln == uniq id lists
		'outall'	: '{0:s}.all'.format(infile),

		#_name of output of dna strings (turned off on 20150417)
		# 'outdna'	: '{0:s}.seq'.format(infile),

		}


#############################################################################_80
#_main_#########################################################################
################################################################################


def run_main(fname=0, outidx='index.txt', **kw):
	import fileinput
	import os
	import itertools

	if os.path.exists(outidx):
		os.unlink(outidx)

	#_open input file
	f = fileinput.input(fname)

	#_keep list of matches and which line they're on in the output file
	dict_out = {}		#_dna_str -> line of output file to store
	dict_idx = {}		#_dna_str -> current collection of uniq_ids

	list_id = []
	dict_sq = {}

	#_loop over two lines of input
	for line0, line1 in itertools.izip_longest(*[f]*2):
		uniq_id = line0.replace('\n','').replace('>','')#_don't .strip()
		dna_str = line1.strip()

		#_build list of unique ids
		idx_id = len(list_id)	#_get current id index
		list_id.append(uniq_id)

		#_build dictionary of
		if dna_str in dict_sq:
			dict_sq[dna_str].append(idx_id)

		else:
			dict_sq[dna_str] = [idx_id]

	#_export two files
	# 1) all unique dna strings (turned off on 20150417)
	# 2) all dna strings and location in original file
	# write_dna(dict_sq, **kw) # (turned off on 20150417)
	write_all(list_id, dict_sq, **kw)

	#_close input
	f.close()


#############################################################################_80
#_end_main_#####################################################################
################################################################################


def write_all(ids, dna, outall='outall.file', sep=' ::: ', **kwargs):
	''' write output file or update it '''
	from os import rename
	import re

	fmt = '{0:s}' + sep + '{1:s}\n'
	with open(outall, 'w') as f:
		for dna_seq, idx in dna.iteritems():
			#_build list of uniq ids
			idx = '; '.join([ids[i] for i in idx])
			f.write(fmt.format(dna_seq, idx))


# (turned off on 20150417)
# def write_dna(strings, outdna='outdna.file', **kwargs):
# 	''' write list of unique strings '''
# 	with open(outdna, 'w') as f:
# 		[f.write('{0:s}\n'.format(l)) for l in strings]


if __name__ == '__main__':
	run_main(**kwargs)
