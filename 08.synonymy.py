#!/usr/local/bin/python3

# Author : Cameron Prybol
# Email : cameron.prybol@gmail.com
# Created : 2014.07.02
# updated :
# description : reports which vcf file entries create non-synonymous encoding sequences

import math
import sys
import subprocess
import difflib

GFF = sys.argv[1]	# gff file
TF = sys.argv[2]	# transcript file
MUT_VCF = sys.argv[3]	# mutant vcf file
OUT_FILE_BASE = sys.argv[4]
DNA_codon_table = {
#         T             C             A             G
# T
        'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT':'C', # TxT
        'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC':'C', # TxC
        'TTA': 'L', 'TCA': 'S', 'TAA': '-', 'TGA':'-', # TxA
        'TTG': 'L', 'TCG': 'S', 'TAG': '-', 'TGG':'W', # TxG
# C
        'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT':'R', # CxT
        'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC':'R', # CxC
        'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA':'R', # CxA
        'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG':'R', # CxG
# A
        'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT':'S', # AxT
        'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC':'S', # AxC
        'ATA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA':'R', # AxA
        'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG':'R', # AxG
# G
        'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT':'G', # GxT
        'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC':'G', # GxC
        'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA':'G', # GxA
        'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG':'G' # GxG
}
def translate_DNA_codon(codon):
        return DNA_codon_table[codon]

revcom_codon_table = {
	'A' : 'T',
	'C' : 'G',
	'G' : 'C',
	'T' : 'A',
	'N' : 'N'	# if N for aNy base, leave as N
	}

def revcom(base):
	return revcom_codon_table[base]


###########################################################################################
#	This section reads through the filtered VCF file and parses the contig, bp-position
#	reference base and alternate base
###########################################################################################

mut_vcf_values = [ [] for i in range(7)]

with open(MUT_VCF) as vcf_file:
	for line in vcf_file:
		line = line.strip("\n")
		if line.startswith("Supercontig"):
			values = line.split()
			contig = int(values[0].split(".")[1])-1		# converted to base 0
			pos = int(values[1])-1				# converted to base 0
			ref = values[3]
			alt = values[4]
			extra_fields = values[7].split(';')
			read_depth = ''
			DP4 = ''
			for field in extra_fields:
				if field.startswith('DP='):
					read_depth = field
				elif field.startswith('DP4='):
					DP4 = field
			mut_vcf_values[contig].append( [ pos , ref , alt , read_depth, DP4] )
	
vcf_file.close()

###########################################################################################
#	Reads through supplied GFF file and extracts feature name, start site, stop site
#	the CDS sequence, if it's on the forward strand, and the reading frame	
###########################################################################################

gene_list = [ [] for i in range(7)]
	
with open(GFF) as f:
	next(f)			#skip over gff header line
	current_gene = ''	# hold ID of current gene
	current_index = ''	# hold ID of current gene
	for line in f:
		line = line.strip("\n")
		values = line.split()

		# establish a 2d array for each gene
		if values[2] == 'gene':
			contig = int(values[0].split(".")[1]) - 1	# use base 0
			start = int(values[3])-1
			stop = int(values[4])-1
			parent = values[8].split(";")[0].split("=")[1]
			current_gene = parent
			gene_list[contig].append([parent,start,stop,[],[]])
			current_index = len(gene_list[contig]) - 1

		else:
			feature_type = values[2]
			contig = int(values[0].split(".")[1]) - 1	# use base 0
			start = int(values[3])-1
			stop = int(values[4])-1
			parent = values[8].split(";")[1].split("=")[1][0:8]
			strand = values[6]
			frame = values[7]
			if gene_list[contig][current_index][0] == parent:
				gene_list[contig][current_index][3].append([feature_type,start,stop,strand,frame])
			else:
				print('error, attempting to write feature to incorrect parent gene')


#############################################################################################
#	run through each gene to find which ones contain VCF snps
#############################################################################################

filtered_gene_list = [ [] for i in range(7)]


for i in range(0,7):
	vcf_contig = mut_vcf_values[i]
	if vcf_contig:	# check if any snps are in that contig, thus skipping contigs without any snps
		gff_contig = gene_list[i]
		#	for each gene, run through the vcf list and check to see if it has a snp located on it
		#	if so, add it to the filtered_gene_list
		for gene in gff_contig:
			start = gene[1]
			stop = gene[2]
			contains_snp = 0
			for snp in vcf_contig:
				snp_start = snp[0]
				snp_len = len(snp[1])
				snp_stop = snp_start + snp_len
				if ((start <= snp_start <= stop) or (start <= snp_stop <= stop)):
					contains_snp = 1
					gene[4].append(snp)
			if contains_snp == 1:
				filtered_gene_list[i].append(gene)

#############################################################################################
#       write list of effected genes to file
#############################################################################################

f_name = OUT_FILE_BASE + '.gene_list'
f = open(f_name,'w')

for contig in filtered_gene_list:
	for gene in contig:
		f.write(gene[0] + '\n')

f.close()



######################################################################################################################
#	This section reads through the transcript file and creates a hash directory for the program to use later
#	in order to create mutant sequences from the original sequence
######################################################################################################################

with open(TF) as tf_file:
	d = []
	CONTIG = ''
	sequence = ''
	for line in tf_file:
		line=line.strip("\n")

		# header line
		if line.startswith(">"):
			###########################################################################
			#	if not the first entry, then CONTIG and sequence should already have
			#	text in them (python evaluates this to True). Check to confirm, and if true, then write
			#	the contents to the dictionary, then reset values
			###########################################################################
			if (CONTIG and sequence):
				d.append([CONTIG,sequence])
				sequence = ''

			# >Supercontig_12.1 of Neurospora crassa OR74A
			CONTIG = line.split()[0][1:]
		# sequence data line
		else:
			sequence += line


#need to write final values to dictionary
d.append([CONTIG,sequence])


###########################################################################################
#	assign each vcf to a sub-feature of the gene
###########################################################################################

snp_info = []

for i in range(0,7):
	contig = filtered_gene_list[i]
	if contig:
		for gene in contig:
			parent = gene[0]
			sub_feature_list = gene[3]
			snp_list = gene[4]
			for snp in snp_list:
				snp_start = snp[0]
				snp_len = len(snp[1])
				snp_stop = snp_start + snp_len
				matching_sub_features = []
				for feature in sub_feature_list:
					start = feature[1]
					stop = feature[2]
					if (((start <= snp_start <= stop) or (start <= snp_stop <= stop)) and feature[0]!='mRNA'):
						matching_sub_features.append(feature)
				if not matching_sub_features:
					matching_sub_features.append('intron')
				info = snp
				info.append(parent)
				info.append('Supercontig=' + str(i+1))
				info.append(matching_sub_features)
				snp_info.append(info)


f_name = OUT_FILE_BASE + '.snp_mapping'
f = open(f_name,'w')

for snp in snp_info:
	features = []
	for feat in snp[7]:
		if ((feat[0] != 'exon') and (feat != 'intron')): 
			features.append(' '.join(list(map(str, feat))))
		if (feat == 'intron'):
			features.append(feat)
	temp = '\t'.join(list(map(str, snp[0:7]))) + '\t' + '\t'.join(list(map(str, features)))
	f.write(temp + '\n')

f.close()

########################################################################################################################
#	translate each CDS region to check for synonymy
########################################################################################################################

CDS_list = []

f_name = OUT_FILE_BASE + '.translations'
f = open(f_name,'w')

# step 1: flip data structure around to have snps listed by CDS, instead of CDSs listed by snp
for snp in snp_info:
	for feature in snp[7]:
		if feature[0] == 'CDS':
			present = 0
			index = ''
			for cds in CDS_list:
				if cds[0] == feature:
					present = 1
					index = CDS_list.index(cds)
			if present == 0:
				CDS_list.append([feature, [] ])
				index = len(CDS_list) - 1
			CDS_list[index][1].append(snp[0:7])

# step 2: obtain sequence for the CDS
for cds in CDS_list:
	contig = int(cds[1][0][6].split('=')[1])-1	# convert to base 0
	start =  cds[0][1]
	stop =  cds[0][2]
	forward = 1
	if cds[0][3] == '-':
		forward = 0
	phase = int(cds[0][4])
	seq = d[contig][1][start:stop+1]	# add +1 to include final base

# step 3: process snps to create mutant DNA, if reverse, reverse complement the DNA, and shift my phase
	mut_seq = ''
	for snp in reversed(cds[1]):
		pos = snp[0]
		off_set = pos - start
		mut_seq = seq[:off_set] + snp[2] + seq[off_set+len(snp[1]):]
	
	###### test to assure correct replacement (it works!) ##############
#		print(seq[off_set-2:off_set+len(snp[1])+2], mut_seq[off_set-2:off_set+len(snp[2])+2], snp[1], snp[2],sep='\t')
	
	if forward == 0:
		seq = seq[::-1]
		temp = ''
		for i in seq:
			temp += revcom(i)
		seq = temp
		mut_seq = mut_seq[::-1]
		temp = ''
		for i in mut_seq:
			temp += revcom(i)
		mut_seq = temp
	seq = seq[phase:]
	mut_seq = mut_seq[phase:]

# step 4: translate and check if AA sequenes are synonymous
	ref_AA = ''
	mut_AA = ''
	for i in range(0,(len(seq)//3)*3,3):
		try:
			ref_AA += translate_DNA_codon(seq[i:i+3])
		except:
			print('failed to transate')
	for i in range(0,(len(mut_seq)//3)*3,3):
		try:
			mut_AA += translate_DNA_codon(mut_seq[i:i+3])
		except:
			print('failed to transate')
	if (ref_AA != mut_AA):
		text = [cds[1][0][5], cds[0][1], cds[0][2], cds[0][3], cds[0][4], ':']
		text = ' '.join(list(map(str,text)))
		f.write(text + '\n')
		for line in (difflib.Differ().compare(ref_AA.splitlines(1),mut_AA.splitlines(1))):
			f.write(line + '\n')
	else:
		text = [cds[1][0][5], cds[0][1], cds[0][2], cds[0][3], cds[0][4], ':', 'synonymous']
		text = ' '.join(list(map(str,text)))
		f.write(text + '\n')
			
f.close()
sys.exit()
