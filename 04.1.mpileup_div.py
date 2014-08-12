#!/usr/local/bin/python3

# Author : Cameron Prybol
# email : cameron.prybol@gmail.com
# date : 2014.07.17
# description: for each contig, takes in 2 mpileup files, and joins them while appending a calculated column with the ratio of the counts for file2/file1

import sys

lens = [ [] for i in range(0,7) ]

f1_pre = sys.argv[1]	# full path for original Oak Ridge mpileup file
f2_pre = sys.argv[2]	# full path for snpless Oak Ridge mpileup file
out_pre = sys.argv[3]	# prefix to name output files, will append .Supercontig_12.(1-7) to end
contig_list = []

for i in range(1,8):
	contig_list.append('Supercontig_12.' + str(i))

##################################################################################################################
#	for each contig, for each position, determine count at position
##################################################################################################################
for i in range(0,len(contig_list)):
	OR_counts = []
	current_pos = 1
	f1_name = f1_pre + "." + contig_list[i]
	with open(f1_name,'r') as f1:
		for line in f1:
			line = line.strip("\n")
			pos = int(line.split()[1])
			while (current_pos < pos):	# fill in positions without values
				OR_counts.append([current_pos,0])
				current_pos += 1
			count = int(line.split()[3])
			OR_counts.append([pos,count])
			current_pos += 1
	f1.close()
	current_pos = 1
##################################################################################################################
#	for each contig, for each position, determine count at position
##################################################################################################################
	f2_name = f2_pre + "." + contig_list[i]
	with open(f2_name,'r') as f2:
		for line in f2:
			line = line.strip("\n")
			pos = int(line.split()[1])
			while (current_pos < pos):	# fill in positions without values
				OR_counts[current_pos-1].append(0)
				current_pos += 1
			count = int(line.split()[3])
			OR_counts[current_pos-1].append(count)
			current_pos += 1
	f2.close()
	out_name = out_pre + "." + contig_list[i]
	out = open(out_name,'w')
##################################################################################################################
#	calculate ratio for counts with snps_removed vs. original, and write all position contents to file
##################################################################################################################
	for n in range(0,len(OR_counts)):
		value = 0
		if (OR_counts[n][1] > 0):
			value = OR_counts[n][2] / OR_counts[n][1]
		OR_counts[n].append(value)
		temp = contig_list[i]+"\t"+str(OR_counts[n][0])+"\t"+str(OR_counts[n][1])+"\t"+str(OR_counts[n][2])+"\t"+str(OR_counts[n][3])+"\n"
		out.write(temp)
	out.close()
