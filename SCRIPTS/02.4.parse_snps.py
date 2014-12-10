#!/usr/local/bin/python3

# Author: Cameron Prybol
# Description: parses SNPS out of mauriceville mapped VCF file

import sys
import subprocess
import pandas as pd

#print(sys.argv[1])
INFILE = sys.argv[1]
OUTFILE = sys.argv[2]
out = open(OUTFILE,'w')

vcf_data = pd.io.parsers.read_csv(INFILE,sep="\t",compression='gzip',header=49)
vcf_data = vcf_data.drop(['FORMAT','/escratch4/cprybol1/cprybol1_Nov_19/MV_MAP_BAM/mv_sim.sorted.bam'],axis=1)
vcf_data['#CHROM'] = vcf_data['#CHROM'].apply(lambda x: str(x).split(".")[1])
print(vcf_data)
