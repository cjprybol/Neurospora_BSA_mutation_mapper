#!/usr/local/bin/python3

# Author : Cameron Prybol
# Email : cameron.prybol@gmail.com
# Created : 2014.07.02
# updated :
# description : reports which vcf file entries create non-synonymous encoding sequences

import math
import sys

GFF = sys.argv[1]	# gff file
SAM = sys.argv[2]	# SAM
FILTER = sys.argv[3]
OUT_FILE_BASE = sys.argv[4]

###########################################################################################
#	Reads through supplied GFF file and extracts feature name, start site, stop site
#	the CDS sequence, if it's on the forward strand, and the reading frame
###########################################################################################

gene_list = [ [] for i in range(7)]

with open(GFF) as f:
	next(f)                 #skip over gff header line
	current_gene = ''       # hold ID of current gene
	current_index = ''      # hold ID of current gene
	for line in f:
		line = line.strip("\n")
		values = line.split()

		# establish a 2d array for each gene
		if values[2] == 'gene':
			contig = int(values[0].split(".")[1]) - 1       # use base 0
			start = int(values[3])-1
			stop = int(values[4])-1
			parent = values[8].split(";")[0].split("=")[1]
			current_gene = parent
			gene_list[contig].append([parent,start,stop,[],[]])
			current_index = len(gene_list[contig]) - 1

		else:
			feature_type = values[2]
			if ((feature_type != 'mRNA') and (feature_type != 'CDS')):
				contig = int(values[0].split(".")[1]) - 1       # use base 0
				start = int(values[3])-1
				stop = int(values[4])-1
				parent = values[8].split(";")[1].split("=")[1][0:8]
				strand = values[6]
				frame = values[7]
				if gene_list[contig][current_index][0] == parent:
					gene_list[contig][current_index][3].append([feature_type,start,stop,strand,frame])
				else:
					print('error, attempting to write feature to incorrect parent gene')

#	for contig in gene_list:
#		for gene in contig:
#			print(gene)


#############################################################################################
#	pull out SA reads from the sam file
#############################################################################################

SA_info = [ [] for i in range(7) ]

with open(SAM) as f:
	for line in f:
		line.strip("\n")
		values = line.split()
		SA = ''
		for value in values:
			if value.startswith('SA:Z:'):
				SA = value[5:-1].split(',')[0:3]
		name = values[0]
		contig1 = int(values[2].split('.')[1])-1	# use 0-base index
		start1 = int(values[3]) - 1		# use 0-base index
		stop1 = start1 + len(values[9])
		read = values[9]
		contig2 = int(SA[0].split('.')[1])-1	# use 0-base index
		if SA[2]=='+':
			start2 = int(SA[1]) - 1		# use 0-base index
			stop2 = start2 + len(values[9])
		elif SA[2]=='-':	# if reverse strand, need to flip start and stop
			stop2 = int(SA[1]) - 1
			start2 = stop2 - len(values[9])
		else:
			print('error')
		# correct for SA's that map to one of the extra contigs
		if (contig1 < 7):
			SA_info[contig1].append([start1,stop1,read,name])
		if (contig2 < 7):
			SA_info[contig2].append([start2,stop2,read,SA[2],name])



#############################################################################################
#	run through each gene to find which ones contain SA reads
#############################################################################################

filtered_gene_list = [ [] for i in range(7)]


for i in range(0,7):
	SA_contig = SA_info[i]
	if SA_contig:	# check if any SAs are in that contig, thus skipping contigs without any SA reads
		gff_contig = gene_list[i]
		#	for each gene, run through the SA list and check to see if it has a SA read located on it
		#	if so, add it to the filtered_gene_list
		for gene in gff_contig:
			start = gene[1]
			stop = gene[2]
			contains_SA = 0
			for sa_read in SA_contig:
				sa_start = sa_read[0]
				sa_stop = sa_read[1]
				if ((start <= sa_start <= stop) or (start <= sa_stop <= stop)):
					contains_SA = 1
					gene[4].append(sa_read)
			if contains_SA == 1:
				filtered_gene_list[i].append(gene)

###########################################################################################
#	assign each sa-read to a sub-feature of the gene
###########################################################################################

SA_mapping = []

for i in range(0,7):
	contig = filtered_gene_list[i]
	if contig:
		for gene in contig:
			parent = gene[0]
			gene_start = gene[1]
			gene_stop = gene[2]
			sub_feature_list = gene[3]
			SA_list = gene[4]
			for SA in SA_list:
				SA_start = SA[0]
				SA_stop = SA[1]
				matching_sub_features = []
				for feature in sub_feature_list:
					start = feature[1]
					stop = feature[2]
					if ((start <= SA_start <= stop) or (start <= SA_stop <= stop)):
						matching_sub_features.append(feature)
				if not matching_sub_features:
					matching_sub_features.append('intron')
#				print(SA)
				info = SA[0:3]
				info.append(SA[-1])
				info.append(parent)
				info.append('Supercontig=' + str(i+1))
				info.append(gene_start)
				info.append(gene_stop)
				info.append(matching_sub_features)
				SA_mapping.append(info)



##############################################################################################
#	import filter list and filter any chimeric genes outside of windows of interest
##############################################################################################

filter_list = open(FILTER,'r')

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


filtered_SA_list = []

#print(len(SA_mapping),'reads in genes')

for sa in SA_mapping:
#	print(sa)
	data = [sa[3],sa[0],sa[4],sa[5]]
	for i in sa[8]:
		if (i != 'intron'):
			temp = ' '.join(map(str,i[0:3]))
		else:
			temp = i
		data.append(temp)
	contig = int(sa[5].split('=')[1]) - 1 # base 0
	start = int(sa[6]/1000)
	stop = int(sa[7]/1000)
	if ((start in SC_locations[contig]) or (stop in SC_locations[contig])):
		filtered_SA_list.append(data)

set_list = []
for SA in filtered_SA_list:
	temp = '\t'.join(map(str,SA))
	set_list.append(temp)

set_list = list(set(set_list))

set_list.sort()
current_read = ''
f_name = OUT_FILE_BASE + '.SA_list'
f = open(f_name,'w')
for set in set_list:
	if (current_read and (current_read != set.split()[0])):
		f.write('\n')
		current_read = set.split()[0]
	else:
		current_read = set.split()[0]
	f.write(set + '\n')

f.close()
