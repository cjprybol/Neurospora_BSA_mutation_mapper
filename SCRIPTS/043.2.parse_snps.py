#!/usr/local/bin/python3

import sys
import subprocess
import pandas as pd

INFILE = sys.argv[1]

OUTFILE = sys.argv[2]

head = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','OTHER']

vcf_data = pd.io.parsers.read_csv(INFILE,sep="\t",header=None)

vcf_data.columns = head

def DP4(info_list):
	temp = str(info_list).split(";")
	return_list = list(filter(lambda item: item.startswith("DP4"), temp))[0]
	return return_list

# pull out DP4 info for easier use
vcf_data['DP4'] = vcf_data['INFO'].apply(lambda x: DP4(x))

def RATIO(DP4_data):
	temp = DP4_data.split("=")[1].split(",")
	ref_total = int(temp[0]) + int(temp[1])
	alt_total = int(temp[2]) + int(temp[3])
	return alt_total/(ref_total + alt_total)

# calculate ratio of reads matching reference (Oak Ridge) vs. matching strain (Mauriceville)
vcf_data['RATIO'] = vcf_data['DP4'].apply(lambda x: RATIO(x))
	
# drop all snps with quality less than 25
vcf_data = vcf_data[(vcf_data['QUAL'] >= 25)]

# drop all snps less than 80%
vcf_data = vcf_data[(vcf_data['RATIO'] > .8)]

def COUNT(DP4_data):
	temp = DP4_data.split("=")[1].split(",")
	ref_total = int(temp[0]) + int(temp[1])
	alt_total = int(temp[2]) + int(temp[3])
	return (ref_total + alt_total)

vcf_data['COUNT'] = vcf_data['DP4'].apply(lambda x: COUNT(x))
	
vcf_data.sort(['#CHROM', 'POS'], inplace=True)
	
vcf_data.reset_index(drop=True, inplace=True)

vcf_data.drop(['DP4', 'RATIO', 'COUNT'], axis=1, inplace=True)

#print(vcf_data)
# save to tab seperated file
vcf_data.to_csv(OUTFILE,sep="\t",index=False, header=False)
