#!/usr/local/bin/python3

# Author : Cameron Prybol
# email : cameron.prybol@gmail.com
# date : 2014.07.24
# description : takes in a .sam file, a list of filter start sites, and an output filename
#	the program reads the filter start site list and removes all reads in the sam file not
#	within the specified regions. regions are 100kb windows
###############################################################################
#	filter site text file format
###############################################################################
#	12.1
#	0-300
#	3300-10800
#	
#	12.3
#	3500-3700
#	
#	12.4
#	...
##############################################################################

import sys
import math

filter_list = open(sys.argv[2],'r')
out_file = open(sys.argv[3],'w')

SC_locations = [[] for i in range(7)]

contig = ''
for line in filter_list:
	line = line.strip("\n")
	if (line.startswith("12.")):
		contig = int(line.split(".")[1])
	elif (len(line) != 0):
		line=line.strip("\n")
		# example: 1970-8117
		start = int(line.split("-")[0])
		end = int(line.split("-")[1])
		for x in range(start,end+1):
			SC_locations[contig-1].append(x)


empty = 1
for supercontig in SC_locations:
	if supercontig:
		empty = 0
	
if (empty == 1):
	print('filter file is empty, aborting this filter')
	sys.exit()

samfile = open(sys.argv[1],'r')
for datum in samfile:
	temp = datum.split()
	kb = 0
	if (len(temp[3]) > 3):
		kb = int(temp[3][:-3])
	try:
		if (int(temp[2].split(".")[1]) <= 7):
			contig = int(temp[2].split(".")[1])
			if (kb in SC_locations[contig-1]):
				out_file.write(datum)
	except:
		print("unexpected .sam line encountered")
		print(datum)

samfile.close()
filter_list.close()
out_file.close()
