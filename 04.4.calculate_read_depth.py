#!/usr/local/bin/python3

import sys

filename = sys.argv[1]
print(filename)

num_lines = 0
total_coverage = 0
with open(filename) as f1:
	for line in f1:
		num_lines += 1
		line = line.strip('\n').split()
		total_coverage += int(line[3])
#		print(line)

print(total_coverage/num_lines)

