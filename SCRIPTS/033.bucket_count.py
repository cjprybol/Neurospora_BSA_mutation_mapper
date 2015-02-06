#!/usr/local/bin/python3

import sys
import pandas as pd
import csv

IN_FILE = sys.argv[1]
OUT_FILE = sys.argv[2]
WINDOW = int(sys.argv[3])
SLIDE = int(sys.argv[4])

window = WINDOW
slide = SLIDE

snp_data = pd.read_table(IN_FILE,sep="\t",header=0)

out_file = open(OUT_FILE,'w')
writer = csv.writer(out_file, delimiter="\t")

writer.writerow( [ 'CONTIG', 'KB', 'TOTAL', 'RATIO', 'MISMATCH_RATIO'] )

# for each chromosome, in order
for i in range (1,8):

# pull out all snps that are on that chromosome
	contig_data = snp_data[snp_data['CONTIG'] == i]
# create a new row containing which KB the snp is in
	contig_data['KB'] = contig_data['POS']//1000

# the list is sorted, so the smallest KB is in the first entry
	min_kb = contig_data.iloc[0]['KB']
# and the largest is in the last entry
	max_kb = contig_data.iloc[len(contig_data.index)-1]['KB']

# for each slide from the lowest KB containing a snp to the highest KB containing a snp
	for window_start in range(min_kb,max_kb+1,slide):
		window_end = window_start + window
	# subset all data by finding the snps that fall between the window start and end
		window_data = contig_data[(contig_data['KB'] >= window_start) & (contig_data['KB'] < window_end)]
	# find the total number of reads matching the reference or alternate SNP in the window
		total = window_data['REF'].sum() + window_data['ALT'].sum()
		ratio = 0
		entries = len(window_data.index)
	# for each snp in the window calculate the normalized ratio, and add it to the total window ratio
	#	normalized based on the number of reads at that location compared to the total reads in the window
	#	thus giving more weight to snp loci with more reads
		for j in range(0,entries):
			ref = window_data.iloc[j]['REF']
			alt = window_data.iloc[j]['ALT']
			entry_total = ref+alt
			if (entry_total > 0):
				entry_ratio = ref / entry_total
				normalized_entry_ratio = entry_ratio * ( entry_total / total )
				ratio += normalized_entry_ratio
		mismatch_total = window_data['MIS'].sum()
		mis_ratio = ''
	# try except to handle snps with no reads (division by 0)
		try:
			mis_ratio = mismatch_total / (total + mismatch_total)
		except:
			mis_ratio = "NA"
		if (ratio == 0):
			ratio = "NA"
		writer.writerow( [ i, window_start, total, ratio, mis_ratio ] )
