#!/usr/local/bin/python3

# Author: Cameron Prybol
# Description: parses SNPS out of mauriceville mapped VCF file

import sys
import subprocess
import pandas as pd

INFILE = sys.argv[1]
OUTFILE = sys.argv[2]

# read infile, which is tab seperated, and in .gz compression format, with the header line on line 49
vcf_data = pd.io.parsers.read_csv(INFILE,sep="\t",compression='gzip',header=49)
	
# drop the 'FORMAT' and '/escratch4/cprybol1/cprybol1_Nov_19/MV_MAP_BAM/mv_sim.sorted.bam' columns
vcf_data = vcf_data.drop([vcf_data.columns.values[8],vcf_data.columns.values[9]],axis=1)

# for each x in #CHROM column, split Chromosome '12.($supercontig)' and replace
#	original value with '$supercontig'
vcf_data['#CHROM'] = vcf_data['#CHROM'].apply(lambda x: int(str(x).split(".")[1]))

# reassign vcf_data to only keep supecontigs from 1-7
#	(the others are mostly repeat regions and don't contain genes of interest)
vcf_data = vcf_data[vcf_data['#CHROM'] <= 7]

# only keep specified columns (drop 'ID', 'QUAL', 'FILTER', 'INFO')
vcf_data = vcf_data[['#CHROM', 'POS', 'REF', 'ALT']]

# rename column headers to drop the '#' infront of '#CHROM'
vcf_data.columns = ['CHROM', 'POS', 'REF' , 'ALT']

# save to tab seperated file
vcf_data.to_csv(OUTFILE,sep="\t")
