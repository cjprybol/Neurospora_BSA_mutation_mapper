#!/usr/local/bin/python3

import sys
import subprocess
import pandas as pd
from Bio.Seq import Seq
import difflib

SNP_GFF = sys.argv[1]
FASTA_SEQ = sys.argv[2]
OUTFILE = sys.argv[3]

snp_data = pd.io.parsers.read_table(SNP_GFF, header=None)

snp_data = snp_data[[8, 3, 4, 6, 7, 10, 12, 13]]

def strip_name(x):
	return (x.split("=")[2])

snp_data.columns = ['name', 'start', 'stop', 'strand', 'offset', 'pos', 'ref', 'alt']

snp_data['name'] = snp_data['name'].apply(lambda x: strip_name(x))

snp_data.set_index('name', inplace=True)

seq_data = pd.io.parsers.read_table(FASTA_SEQ, index_col=0, header=0)
seq_data.index.names = ['name']

merged = pd.merge(snp_data, seq_data, left_index=True, right_index=True)

f = open(OUTFILE, 'w')

for i in merged.index:
	position = merged.at[i,'pos'] - merged.at[i,'start']
	snp_len = len(merged.at[i,'ref'])
	offset = merged.at[i,'offset']
	base = ""
	seq = ""
	if (merged.at[i,'strand'] == "+"):
		seq = Seq(merged.at[i,'sequence'])
	elif (merged.at[i,'strand'] == "-"):
		seq = Seq(merged.at[i,'sequence']).reverse_complement()
	new_seq = seq[:position] + merged.at[i,'alt'] + seq[position+snp_len:]
	if (merged.at[i,'strand'] == "-"):
		seq = seq.reverse_complement()
		new_seq = new_seq.reverse_complement()
	ref_tr_seq = str(seq.translate())[offset:]
	new_tr_seq = str(new_seq.translate())[offset:]

	# split every n  nucleotides
	n = 50
	ref_list = [(str(i+1) + "." + ref_tr_seq[i:i+n]) for i in range(0, len(ref_tr_seq), n)]
	new_list = [(str(i+1) + "." + new_tr_seq[i:i+n]) for i in range(0, len(new_tr_seq), n)]

	diff = difflib.Differ().compare(ref_list, new_list)
	expected_AA = (position//3) + 1
	print('Transcript', 'Pos', 'Ref', 'Alt', 'AA Location', sep="\t\t", file=f)
	print(i, merged.at[i,'pos'], merged.at[i,'ref'], merged.at[i,'alt'], expected_AA, sep="\t\t", file=f)
	print(file=f)
	print('sequence:', "\n".join(diff), sep="\n", file=f)
	print(file=f)
