#!/usr/local/bin/python3

# Author: Cameron Prybol
# Created: 2014.06.22
# Last Updated:
# Description: removes duplicate reads from Oak Ridge .sam file
# run in command line as:
#	> python3 9.sub_script.py "full to sorted_with_duplicates file" "full path to outfile"

import sys

INFILE = sys.argv[1]
OUTFILE = sys.argv[2]
out = open(OUTFILE,'w')

read_ID = ""
duplicate = 0
snp_count = 0
exact_match_count = 0
with open(INFILE) as f:
	for line in f:
		temp = line.split()
		if len(temp) == 1:
			read_ID = temp[0]
		elif ((temp[0] == read_ID) and ("NM:i:0" not in line)):
			snp_count += 1
		elif ((temp[0] == read_ID) and ("NM:i:0" in line)):
			exact_match_count += 1
			out.write(line)
		else :
			out.write(line)

print(str(snp_count) + " MV snps removed")
print(str(exact_match_count) + " exact matches untouched")
f.close()
out.close()
