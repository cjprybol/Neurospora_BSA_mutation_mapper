#!/usr/local/bin/python3

# Author : Cameron Prybol
# Email : cameron.prybol@gmail.com
# Created : 2014.06.24
# updated :
# description : creates a bam file only including reads in regions listed in .gff file

import sys
import subprocess

GFF = sys.argv[1]
SAM = sys.argv[2]
OUTFILE = sys.argv[3]

gff_range_set = [ [] for i in range(0,20) ] # instantiate 20 lists, one for each supercontig


###########################################################################################
#	populate a list of lists, one per supercontig, full of base pairs within a gff feature
###########################################################################################

with open(GFF) as f:
	next(f)		#skip over gff header line
	for line in f:
		temp = line.split()
		gff_feature_supercontig = int(temp[0].split(".")[1]) - 1
		gff_feature_start = int(temp[3])
		gff_feature_end = int(temp[4])
		for i in range( gff_feature_start , gff_feature_end + 1 ): 
			gff_range_set[gff_feature_supercontig].append(i)
f.close()

###########################################################################################
#	for each supercontig range set, keep only unique values
###########################################################################################

for i in range(0,len(gff_range_set)):
	if gff_range_set[i]:	#check to see if non-empty
		gff_range_set[i] = set(gff_range_set[i])

out = open(OUTFILE,'w')

###########################################################################################
#	for each read, check to see if that read overlaps with a gff feature
#	if so, write that read and it's full .sam file line to a new .sam file
###########################################################################################

counter = 0
current = 0
with open(SAM) as f:
	for line in f:
		counter = counter + 1
		temp = line.split()
		contig = int(temp[2].split(".")[1]) - 1
		read_start_pos = int(temp[3])
		gff_feature_flag = 0
		for i in range(read_start_pos, read_start_pos + 52):	#read length = 51, plus 1 because range is final value exclusive
			if i in gff_range_set[contig]:
				gff_feature_flag = 1
				break
		if (gff_feature_flag == 1):
			out.write(line)
f.close()
out.close()


sys.exit()
