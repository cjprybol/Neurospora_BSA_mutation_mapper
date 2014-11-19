#!/usr/local/bin/python3

# Author: Cameron Prybol
# Created: 2014.09.19
# Last Updated:
# Description: removes duplicate snps from vcf file
# run in command line as:
#	> python3 'script_name' "full to sorted_with_duplicates file" "full path to outfile"

import sys

INFILE = sys.argv[1]
OUTFILE = sys.argv[2]
out = open(OUTFILE,'w')


current_uniq = []
uniq_count = 0
with open(INFILE) as f:
	for line in f:
		temp = line.split()
		if len(temp) == 4:
			current_uniq = [temp[0],temp[1]]
		elif current_uniq:
			if ((temp[0] == current_uniq[0]) and (temp[1] == current_uniq[1])):
				uniq_count += 1
				out.write(line)

print(uniq_count)

f.close()
out.close()
