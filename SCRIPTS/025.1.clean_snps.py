#!/usr/local/bin/python3

import sys

INFILE = sys.argv[1]
OUTFILE = sys.argv[2]

#print(INFILE)
#print(OUTFILE)

out = open(OUTFILE,'w')

prior = ''
skipped = 0
with open(INFILE) as f:
	prior = f.readline().split()
	for line in f:
		line = line.split()
		if len(prior) == 5 and prior == line[0:5]:
			skipped += 1
#			print("prior",prior,sep="\t")
#			print("current",line,sep="\t")
#			print(skipped)
		elif len(line) > 5:
			out.write("\t".join(line)+"\n")
		prior = line
#print(skipped)
