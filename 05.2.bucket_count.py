#!/usr/local/bin/python3

# Author : Cameron Prybol
# email : cameron.prybol@gmail.com
# date : 2014.07.17
# description: creates bucketed counts for mpileup counts to allow for graphing in R

import sys


count = 5000	# size of buckets
slide_size = count // 5
in_file = ''
out_file = ''
if (len(sys.argv) == 2):
	in_file = sys.argv[1]
elif (len(sys.argv) == 3):
	in_file = sys.argv[1]
	count = int(sys.argv[2])

out_file = in_file + '_' + str(count) + 'bp_count'

counts = []

for i in range(0,5):
	current_bucket = 0
	offset = i * slide_size 
	current_total = 0
	records_read = 0
	records_skipped = 0
	with open(in_file) as f:
		for q in range(0,offset):
			next(f)
		for line in f:
			line = line.strip("\n").split()
			temp = int(line[1])-1
			bucket = (temp-offset)//count
			if (bucket > current_bucket):
				if (records_read > 0):
					counts.append([(current_bucket * count) + 1 + offset,current_total/records_read,records_read,records_skipped])
				else:
					counts.append([(current_bucket * count) + 1 + offset,'0',records_read,records_skipped])
				current_bucket = bucket
				current_total = 0
				records_read = 0
				records_skipped = 0
			if (int(line[2]) != 0 and int(line[3]) != 0 and float(line[4]) != 0):	#don't pull down average by including positions without any reads
				current_total += float(line[4])
				records_read += 1
			else:
				records_skipped += 1
	if (records_read > 0):
		counts.append([(current_bucket * count) + 1 + offset,current_total/records_read,records_read,records_skipped])
	else:
		counts.append([(current_bucket * count) + 1 + offset,'0',records_read,records_skipped])

counts.sort()

out = open(out_file,'w')
	
for value in counts:
	output = '\t'.join(map(str, value)) + '\n'
	print(output.strip())
	out.write(output)

out.close()