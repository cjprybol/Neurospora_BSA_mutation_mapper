#!/usr/local/bin/python3

# Author : Cameron Prybol
# Email : cameron.prybol@gmail.com
# Created : 2014.07.16
# updated :
# description : filter mpileup to only retain positions within GFF features

import sys

GFF = sys.argv[1]
prefix = sys.argv[2]
#OUTFILE = sys.argv[3]

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


###########################################################################################
#	if mpileup position is in a GFF feature, then write to file
###########################################################################################

for i in range(1,8):
	in_file = prefix + '.' + str(i)
	out_file = in_file + '.gff_filtered'
	out = open(out_file,'w')
	print("filtering supercontig",i)
	with open(in_file) as f:
		for line in f:
			values = line.strip("\n").split()
			if (int(values[1]) in gff_range_set[i-1]):
				out.write(line)
	f.close()
	out.close()


sys.exit()
