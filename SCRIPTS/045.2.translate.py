#!/usr/local/bin/python3

import sys
import subprocess
import pandas as pd

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

print(merged)
