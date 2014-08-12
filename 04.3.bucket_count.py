#!/usr/local/bin/python3

# Author : Cameron Prybol
# email : cameron.prybol@gmail.com
# date : 2014.07.17
# description: creates bucketed counts for mpileup counts to allow for graphing in R

import sys


count = 5000	# size of buckets
in_file = ''
out_file = ''
if (len(sys.argv) == 2):
	in_file = sys.argv[1]
elif (len(sys.argv) == 3):
	in_file = sys.argv[1]
	count = int(sys.argv[2])

out_file = in_file + '_' + str(count) + 'bp_count'

counts = []

current_bucket = 0
current_total = 0
records_read = 0
with open(in_file) as f:
	for line in f:
		line = line.strip("\n").split()
		temp = int(line[1])-1
		bucket = temp//count
		if (bucket > current_bucket):
			counts.append([(current_bucket * count) + 1,current_total/records_read])
			current_bucket = bucket
			current_total = 0
			records_read = 0
		current_total += float(line[4])
		records_read += 1
counts.append([(current_bucket * count) + 1,current_total/records_read])

out = open(out_file,'w')


for value in counts:
	output = str(value[0]) + "\t" + str(value[1]) + "\n"
	out.write(output)

out.close()
