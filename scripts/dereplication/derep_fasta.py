#!/usr/bin/env python

'''
Authors: Walter Sessions, Jimmy O'Donnell
Find and remove duplicate DNA sequences from a fasta file
usage: python ./this_script.py infile.fasta 'ID1_' derep.fasta derep.map
'''

import os
from collections import Counter
import fileinput
import itertools
import argparse

parser = argparse.ArgumentParser(
	description = 'Remove and count duplicate sequences in a fasta file', 
	)

parser.add_argument('-i', '--fasta_in',
	help = 'File path for input file. Must be a fasta file with *no* wrapped lines!', 
	required = True)

parser.add_argument('-s', '--sample_prefix',
	help = 'First part of string identifying samples (e.g. "ID1=").', 
	required = True)

parser.add_argument('-f', '--fasta_out',
	help = 'File path for output fasta.', 
	required = True)

parser.add_argument('-m', '--mapfile',
	help = 'File path for output mapfile.', 
	required = True)

args = vars(parser.parse_args())

infile = args['fasta_in']
sample_id_start = args['sample_prefix']
outfasta = args['fasta_out']
outmap = args['mapfile']

#############################################################################_80
#_main_#########################################################################
################################################################################


def run_main(fname, sampleID, out_f, out_m):

	#_open input file
	f = fileinput.input(fname)

	list_id = []
	dict_uniqseq = {}

	# sample_pattern = re.compile(samplestringID + "(.*)", re.flags)

	#_loop over two lines of input
	for line0, line1 in itertools.izip_longest(f,f):
		seq_id = line0.replace('\n','').replace('>','')#_don't .strip()
		sample_id = sampleID + '{0:s}'.format(seq_id.split(sampleID)[1])
		dna_str = line1.strip()

		#_build list of unique ids
		idx_id = len(list_id)	#_get current id index
		list_id.append(sample_id) # seq_id

		#_build dictionary of
		if dna_str in dict_uniqseq:
			dict_uniqseq[dna_str].append(idx_id)

		else:
			dict_uniqseq[dna_str] = [idx_id]

	keys_by_length = sorted(dict_uniqseq, 
						key=lambda k: len(dict_uniqseq[k]), reverse = True)
	write_fasta(dict_uniqseq, keys_by_length, out_f)
	write_map(list_id, dict_uniqseq, keys_by_length, out_m)

	#_close input
	f.close()


#############################################################################_80
#_end_main_#####################################################################
################################################################################

def write_fasta(dna_dict, sorted_keys, fasta_output = 'fasta_output.file'):
	''' write fasta file of unique sequences '''
	with open(fasta_output, 'w') as f:
		for index, key in enumerate(sorted_keys):
			# if len(dna_dict[key]) > 1:
			f.write('>DUP_{0:n}'.format(index+1) +
			        ';size={0:n}\n'.format(len(dna_dict[key])) +
					'{0:s}\n'.format(key))

def write_map(id_list, dna_dict, sorted_keys, map_output = 'map_output.file'):
	'''write two column file of sequence name from input and sequence name in output'''
	with open(map_output, 'w') as f:
		for index, key in enumerate(sorted_keys):
			count_per_sample = Counter([id_list[i] for i in dna_dict[key]]).most_common()
		    	f.write('\n'.join([
					'DUP_{0:n}'.format(index+1) + '\t' +
			        '{0:s}'.format(k) + '\t' +
					'{0:n}'.format(v)
					for k, v in count_per_sample]) + '\n')

if __name__ == '__main__':
	run_main(fname = infile, sampleID = sample_id_start, out_f = outfasta, out_m = outmap)
