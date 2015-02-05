#!/usr/local/bin/python3

# Author: Cameron Prybol
# Description: parses SNPS out of mauriceville mapped VCF file

import sys
import subprocess
import pandas as pd

INFILE = sys.argv[1]

#print(INFILE)
OUTFILE = str(INFILE.split(".")[0]) + ".parsed_snps.out"
#print(OUTFILE)

head = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','JUNK']

vcf_data = pd.io.parsers.read_csv(INFILE,sep="\t",header=None)


vcf_data.columns = head

# drop last 2 columns
vcf_data = vcf_data.drop([vcf_data.columns.values[8],vcf_data.columns.values[9]],axis=1)

# for each x in #CHROM column, split Chromosome '12.($supercontig)' and replace
#	original value with '$supercontig'
vcf_data['#CHROM'] = vcf_data['#CHROM'].apply(lambda x: int(str(x).split(".")[1]))

def DP4(info_list):
	temp = str(info_list).split(";")
	return_list = list(filter(lambda item: item.startswith("DP4"), temp))[0]
	return return_list

# keep only the DP4 info from the additional infomation
vcf_data['INFO'] = vcf_data['INFO'].apply(lambda x: DP4(x))

# only keep specified columns (drop 'ID', 'FILTER')
vcf_data = vcf_data[['#CHROM', 'POS', 'REF', 'ALT','QUAL','INFO']]

# rename column headers to drop the '#' infront of '#CHROM'
vcf_data.columns = ['CONTIG', 'POS', 'REF' , 'ALT','QUAL','INFO']

# drop all snps with quality less than 25
vcf_data = vcf_data[(vcf_data['QUAL'] >= 25)]

def RATIO(DP4_data):
	temp = DP4_data.split("=")[1].split(",")
	ref_total = int(temp[0]) + int(temp[1])
	alt_total = int(temp[2]) + int(temp[3])
	return ref_total/(ref_total + alt_total)

# calculate ratio of reads matching reference (Oak Ridge) vs. matching strain (Mauriceville)
vcf_data['RATIO'] = vcf_data['INFO'].apply(lambda x: RATIO(x))

def COUNT(DP4_data):
	temp = DP4_data.split("=")[1].split(",")
	ref_total = int(temp[0]) + int(temp[1])
	alt_total = int(temp[2]) + int(temp[3])
	return (ref_total + alt_total)

vcf_data['COUNT'] = vcf_data['INFO'].apply(lambda x: COUNT(x))

# drop all locations without perfect consensus
vcf_data = vcf_data[(vcf_data['RATIO'] == 0)]

vcf_data.sort(['CONTIG', 'POS'], inplace=True)
	
# save to tab seperated file
vcf_data.reset_index(drop=True, inplace=True)
#print(vcf_data)
vcf_data.to_csv(OUTFILE,sep="\t",index=False)
