#!/usr/local/bin/python3

import sys

INFILE = sys.argv[1]
OUTFILE = sys.argv[2]

out = open(OUTFILE,'w')

prior = ''
skipped = 0
with open(INFILE) as f:
	prior = f.readline().split()
	for line in f:
		line = line.split()
		# if the previous line was only 5 entries long (length of inserted snp markers that should be removed)
		#	and the first 5 entries of the current snp match that inserted marker,
		#	skip it, thus removing both from the outfile
		if len(prior) == 5 and prior == line[0:5]:
			skipped += 1
		# if the line is a full entry (greater than 5 columns, indicating it was not from the snp marker list to
		#	be removed), and did not match the prior criteria, it is a snp we want to keep, so
		#	write it out to the file
		elif len(line) > 5:
			out.write("\t".join(line)+"\n")
		prior = line
