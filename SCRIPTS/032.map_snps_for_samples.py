#!/usr/local/bin/python3

import pandas as pd
import sys
import csv

IN = sys.argv[1]
OUT_FILE = sys.argv[2]
SNP_FILE = sys.argv[3]


# set column names for pre-cleaned sam data

col_names=['QNAME','CONTIG','POS','SEQ']

# read in data
sam_data = pd.read_csv(IN, delim_whitespace=True,names=col_names)
#sam_data = pd.read_csv(IN, delim_whitespace=True,names=col_names,nrows=139429)

# calculate end position of each read

sam_data['END_POS'] = sam_data.apply(lambda row: row['POS'] + len(row['SEQ']) - 1,axis=1)

snp_data = pd.read_csv(SNP_FILE, delim_whitespace=True,header=0)
#snp_data = pd.read_csv(SNP_FILE, delim_whitespace=True,header=0,nrows=78100)

f = open(OUT_FILE, 'w')
writer = csv.writer(f,delimiter='\t')
print('saving to:', OUT_FILE, sep='\t')

# iterate over each chromosome in turn
for i in range(1,8):

	sam_chrom_sub = sam_data[sam_data.CONTIG == i]
	sam_data = sam_data[sam_data.CONTIG != i]

	snp_chrom_sub = snp_data[snp_data.CONTIG == i]
	snp_data = snp_data[snp_data.CONTIG != i]

	for snp_index, snp in snp_chrom_sub.iterrows():

		# snp_index returns the index of the original dataframe, not the subset. need to figure that out
#		print('snp',snp_index,'of',len(snp_chrom_sub.index))

		ref_match_count = 0
		alt_match_count = 0
		read_mismatch_count = 0

		# try to account for >1 length vcf SNPS
		ref_snp_extra = len(snp['REF']) - 1
		alt_snp_extra = len(snp['ALT']) - 1
		extra_length = ref_snp_extra
		larger_snp = 'REF'
		# choose bigger of the two
		if alt_snp_extra > ref_snp_extra:
			extra_length = alt_snp_extra
			larger_snp = 'ALT'

		sam_snp_sub = sam_chrom_sub[(sam_chrom_sub.POS <= snp['POS']) & (sam_chrom_sub.END_POS >= (snp['POS'] + extra_length))]

		for matched_read_index, matched_read in sam_snp_sub.iterrows():
			offset = snp['POS'] - matched_read['POS']
			base = matched_read['SEQ'][offset:offset+extra_length+1]

			ref_match = 0
			alt_match = 0
			no_match = 0
			# if just a single base, make a simple comparison
			if (extra_length == 0):
				if (base == snp['REF']):
					ref_match = 1
				if (base == snp['ALT']):
					alt_match = 1
				if (ref_match == 0 and alt_match == 0):
					no_match = 1
				if (ref_match == 1 and alt_match == 1):
					writer.writerow(['ERROR!!!!'])
			# if an indel change, it gets a little complex
			if extra_length > 0:
				# check to see if parsed snp matches over size of each SNP
				if (base[0:len(snp['REF'])] == snp['REF']):
					ref_match = 1
				if (base[0:len(snp['ALT'])] == snp['ALT']):
					alt_match = 1
				# if they both match, it must be the larger, so toggle the smaller to not matching
				if (ref_match == 1 and alt_match == 1):
					if (base[0:len(snp[larger_snp])] == snp[larger_snp]):
						if (larger_snp == 'REF'):
							alt_match = 0
						elif (larger_snp == 'ALT'):
							ref_match = 0
						else:
							writer.writerow(['ERROR!!!!'])
				if (ref_match == 0 and alt_match == 0):
					no_match = 1
			if (no_match == 1):
				writer.writerow([ 'base:', base, 'ref:', snp['REF'], 'ref_match:', ref_match, 'alt:', snp['ALT'], 'alt_match:', alt_match, 'NO MATCH' ])
			else:
				writer.writerow([ 'base:', base, 'ref:', snp['REF'], 'ref_match:', ref_match, 'alt:', snp['ALT'], 'alt_match:', alt_match ])
			ref_match_count += ref_match
			alt_match_count += alt_match
			read_mismatch_count += no_match
		writer.writerow(['CONTIG', i ,'POS', snp['POS'], 'ref_match_count:',ref_match_count,'alt_match_count:',alt_match_count,'read_mismatch_count:',read_mismatch_count ])
#		print('CONTIG', i ,'POS', snp['POS'], 'ref_match_count:',ref_match_count,'alt_match_count:',alt_match_count,'read_mismatch_count:',read_mismatch_count,sep="\t" )
