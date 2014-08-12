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
DIM5_VCF = sys.argv[4]	# DIM5 vcf file
OUT_FILE = sys.argv[5]
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

######################################################################################################################
#	This section reads through the transcript file and creates a hash directory for the program to use later
#	in order to create mutant sequences from the original sequence
######################################################################################################################

print("extracting sequences from supercontigs.fasta file...")
with open(TF) as tf_file:
	d = {}	# create an empty dictionary
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
				d[CONTIG] = sequence
				CONTIG = ''
				sequence = ''

			# >Supercontig_12.1 of Neurospora crassa OR74A
			CONTIG = str(int(line.split()[0].split(".")[1])-1)	#make it zero based index
		# sequence data line
		else:
			sequence += line


#need to write final values to dictionary
d[CONTIG] = sequence

###########################################################################################
#	Reads through supplied GFF file and extracts feature name, start site, stop site
#	the CDS sequence, if it's on the forward strand, and the reading frame	
###########################################################################################

print("extracting features from transcripts.gff file...")
gff_CDS = [ [] for i in range(7)]
gff_transcript = [ {} for i in range(7)]
	
with open(GFF) as f:
	next(f)		#skip over gff header line
	for line in f:
		line = line.strip("\n")
		values = line.split()
		
		# establish a dictionary for each mRNA transcript
		if values[2] == "mRNA":
			# ID=[NCU11405T0];Parent=NCU11405
			parent = values[8].split(";")[0].split("=")[1]
			contig = int(values[0].split(".")[1]) - 1
			start = int(values[3])-1
			stop = int(values[4])-1
			gff_transcript[contig][parent] = [[start,stop]]



		# only CDS regions have reading frame information, so extract all information to compare for synonymy
		elif values[2] == "CDS":

			#Supercontig_12.[1-7]
			contig = int(values[0].split(".")[1]) - 1
			start = int(values[3])-1
			stop = int(values[4])-1
			forward = 1
			if (values[6] == '-'):
				forward = 0
			frame = values[7]
			# ID=CDS:NCU11405T0:1;Parent=[NCU11405T0]
			parent = values[8].split(";")[1].split("=")[1]
			# ID=CDS:[NCU11405T0:1];Parent=NCU11405T0
			name = values[8].split(";")[0].split("=")[1][4:]
			full_contig = d[str(contig)]
			CDS_seq = full_contig[start:stop+1]	# add 1 to include last character

			if CDS_seq:
				gff_CDS[contig].append( [ name , start , stop , CDS_seq , forward, frame] ) #[transcript_name] [start] [stop] [sequence] [strand]
			else:
				print("error",CDS_seq)

		# for all other non-gene regions, just get start and stop to see if any SNPs are in them
		elif values[2] != "gene":
			# ID=[NCU11405T0];Parent=NCU11405
			parent = values[8].split(";")[1].split("=")[1]
			name = values[8].split(";")[0].split("=")[1]
			contig = int(values[0].split(".")[1]) - 1
			start = int(values[3])-1
			stop = int(values[4])-1
			gff_transcript[contig][parent].append([name,start,stop])


f.close()


###########################################################################################
#	This section reads through the VCF file for both the mutant and DIM5 strains
#	and creates an organized list to read through later on
###########################################################################################

print("extracting snp data from mutant strain...")
mut_vcf_values = [ [] for i in range(7)]

with open(MUT_VCF) as vcf_file:
	for line in vcf_file:
		# data line
		if line.startswith("Supercontig"):
			line = line.strip("\n")
			values = line.split()
			contig = int(values[0].split(".")[1])-1
			pos = int(values[1])-1
			ref = values[3]
			alt = values[4].split(",")
			if 'X' in alt: alt.remove('X')	#remove the X base added by samtools mpileup
			mut_vcf_values[contig].append( [ pos , ref , alt[0] ] )

vcf_file.close()

print("extracting snp data from dim5 strain...")
dim5_vcf_values = [ [] for i in range(7)]

with open(DIM5_VCF) as vcf_file:
	for line in vcf_file:
		# data line
		if line.startswith("Supercontig"):
			line = line.strip("\n")
			values = line.split()
			contig = int(values[0].split(".")[1])-1
			pos = int(values[1])-1
			ref = values[3]
			alt = values[4].split(",")
			if 'X' in alt: alt.remove('X')	#remove the X base added by samtools mpileup
			dim5_vcf_values[contig].append( [ pos , ref , alt[0] ] )


vcf_file.close()

###########################################################################################
#	assign SNPS to CDS regions to in order to evaluate synonymy across entire CDS region
###########################################################################################

cross_check_CDS = [ [] for i in range(7)]

	
print("matching SNPs to CDS regions...")
for i in range(0,7):
	current_contig = gff_CDS[i]
	for CDS in current_contig:
		#	gff_CDS[contig][#] = [ name , start , stop , CDS_seq , forward]
		mut_SNPs = []
		dim5_SNPs = []
		for snp in mut_vcf_values[i]:
			# if starts within the CDS, and ends within the CDS...
			if ((CDS[1] <= snp[0] <= CDS[2]) and (CDS[1] <= (snp[0]+len(snp[1])) <= CDS[2])):
				mut_SNPs.append(snp)
		for snp in dim5_vcf_values[i]:
			# if starts within the CDS, and ends within the CDS...
			if ((CDS[1] <= snp[0] <= CDS[2]) and (CDS[1] <= (snp[0]+len(snp[1])) <= CDS[2])):
				dim5_SNPs.append(snp)
		if (len(mut_SNPs) != 0 or len(dim5_SNPs) != 0):	# don't waste time reading through SNP-less transcripts later on
			cross_check_CDS[i].append([CDS,mut_SNPs,dim5_SNPs])
CDS_count = 0
for i in range(0,7):
	CDS_count += len(cross_check_CDS[i])
print("CDS_count :",CDS_count)

#################################################
#	feature not currently being used
#################################################
#	cross_check_features = [ [] for i in range(7)]
#	
#	print("matching SNPs to other gff features...")	
#	for i in range(0,7):
#		current_contig = gff_transcript[i]
#		full_contig_set = list(current_contig.values())
#		for parent_transcript in full_contig_set:
#			for j in range(1,len(parent_transcript)):
#				FEATURE = parent_transcript[j]
#				mut_SNPs = []
#				dim5_SNPs = []
#				for snp in mut_vcf_values[i]:
#					if (FEATURE[1] <= snp[0] <= FEATURE[2]):
#						mut_SNPs.append(snp)
#				for snp in dim5_vcf_values[i]:
#					if (FEATURE[1] <= snp[0] <= FEATURE[2]):
#						dim5_SNPs.append(snp)
#				if (len(mut_SNPs) != 0 or len(dim5_SNPs) != 0):	# don't waste time reading through SNP-less transcripts later on
#					cross_check_features[i].append([FEATURE,mut_SNPs,dim5_SNPs])
#	
#	feature_count = 0
#	for i in range(0,7):
#		feature_count += len(cross_check_features[i])
#	print("feature_count :",feature_count)
#################################################
#	end feature
#################################################




###########################################################################################
#	sub mutant SNPs for original sequence, and check for synonymy
###########################################################################################

#for contig in cross_check_CDS:
print("translating DNA sequences to check for snyonymy...")
#keep track of how many times mpileup reference DOES NOT match user determined reference
conflicts = 0
#keep track of how many times mpileup reference DOES match user determined reference
matches = 0
non_syn_snps = [ [ [] , [] , [] , [] ] for i in range(0,7)]

#for each contig 1-7 (offset by 1 to 0-6)
for n in range(0,7):
	#assign that contig to a variable, CDS_contig
	CDS_contig = cross_check_CDS[n]
	# for each GFF CDS feature on that contig
	for FEATURE in CDS_contig:
		mutant_synonymy = 'T'	# assume all SNPs are synonymous
		dim5_synonymy = 'T'	# assume all SNPs are synonymous
		#create a list in the form [feature and all it's details , mutant non-syn snps , dim5 non-syn snps]
		non_syn_feats = [ FEATURE[0] , [] , [] ]
		# take reference sequence from feature details
		reference_seq = FEATURE[0][3]
		# take reading frame from reference details
		frame = int(FEATURE[0][5])


		####################################################################################
		#	mutant
		####################################################################################
		# create variable to hold the mutant sequence as we build it
		mutant_seq = reference_seq
		# variable to hold the last stop position (start at 0)
		prior_stop = 0
		# figure out the number of mutant SNPs in this feature to evaluate
		snp_count = len(FEATURE[1])
		# variable to hold sequence AFTER point mutation
		post = ''
		# variable to keep track of length change due to indels
		offset = 0
		# for each snp in the feature set...
		for q in range(0,snp_count): 
			# assign the snp and it's details to a variable
			snp = FEATURE[1][q]
			# feature start site
			start = FEATURE[0][1]
			# figure out position of snp relative to CDS start site
			pos = snp[0] - start
			# assign the reference base to a variable
			ref_snp = snp[1]
			# assign mpileup determined mutant snp to variable
			mutant_snp = snp[2]


			# TEST#
			if (ref_snp != mutant_seq[pos+offset:pos+offset+len(ref_snp)]):
				conflicts += 1
			else:
				matches+=1
			# END TEST #


			# determine the sequence up until the point of mutation
			pre = mutant_seq[0:pos+offset]			#draws from mut_seq to keep changes
			# determine all sequence that should occur after the mutation
			# making sure to draw from the reference sequence
			post = reference_seq[pos+len(ref_snp):]
			#create mutant sequence from parts
			mutant_seq = pre + mutant_snp + post
			# keep track of length change due to indels
			offset += len(mutant_snp) - len(ref_snp)
			
		# if 'forward' if 0 (false), then need to make revcom of sequence because its on
		#	the reverse strand
		if FEATURE[0][4] == 0:
			temp = reference_seq[::-1]
			reference_seq = ''
			for base in temp:
				reference_seq += revcom(base)
			temp = mutant_seq[::-1]
			mutant_seq = ''
			for base in temp:
				mutant_seq += revcom(base)
			
		# translate reference and mutant AA sequences to see if synonymous
		reference_AA = ''
		mutant_AA = ''
		for i in range(frame,int(len(reference_seq)/3)+1,3):
			reference_AA += translate_DNA_codon(reference_seq[i:i+3])
		for i in range(frame,int(len(mutant_seq)/3)+1,3):
			try:
				mutant_AA += translate_DNA_codon(mutant_seq[i:i+3])
			except:
				print("fail",mutant_seq[i:i+3])
		if reference_AA != mutant_AA:
			mutant_synonymy = 0
			snp_list = FEATURE[1]
			non_synonymous_snp_set = [snp_list,reference_AA,mutant_AA,mutant_seq]
			non_syn_feats[1] = non_synonymous_snp_set




		####################################################################################
		#	DIM 5
		####################################################################################
		# create variable to hold the DIM_5 sequence as we build it
		DIM_5_seq = reference_seq
		# variable to hold the last stop position (start at 0)
		prior_stop = 0
		# figure out the number of DIM_5 SNPs in this feature to evaluate
		snp_count = len(FEATURE[2])
		# variable to hold sequence AFTER point mutation
		post = ''
		# variable to keep track of length change due to indels
		offset = 0
		# for each snp in the feature set...
		for q in range(0,snp_count): 
			# assign the snp and it's details to a variable
			snp = FEATURE[2][q]
			# feature start site
			start = FEATURE[0][1]
			# figure out position of snp relative to CDS start site
			pos = snp[0] - start
			# assign the reference base to a variable
			ref_snp = snp[1]
			# assign mpileup determined DIM_5 snp to variable
			DIM_5_snp = snp[2]


			# TEST#
			if (ref_snp != DIM_5_seq[pos+offset:pos+offset+len(ref_snp)]):
				conflicts += 1
			else:
				matches+=1
			# END TEST #


			# determine the sequence up until the point of mutation
			pre = DIM_5_seq[0:pos+offset]			#draws from mut_seq to keep changes
			# determine all sequence that should occur after the mutation
			# making sure to draw from the reference sequence
			post = reference_seq[pos+len(ref_snp):]
			#create DIM_5 sequence from parts
			DIM_5_seq = pre + DIM_5_snp + post
			# keep track of length change due to indels
			offset += len(DIM_5_snp) - len(ref_snp)
			

		# if 'forward' if 0 (false), then need to make revcom of sequence because its on
		#	the reverse strand
		if FEATURE[0][4] == 0:
			temp = reference_seq[::-1]
			reference_seq = ''
			for base in temp:
				reference_seq += revcom(base)
			temp = DIM_5_seq[::-1]
			DIM_5_seq = ''
			for base in temp:
				DIM_5_seq += revcom(base)
			
		# translate reference and DIM_5 AA sequences to see if synonymous
		reference_AA = ''
		DIM_5_AA = ''
		for i in range(frame,int(len(reference_seq)/3)+1,3):
			reference_AA += translate_DNA_codon(reference_seq[i:i+3])
		for i in range(frame,int(len(DIM_5_seq)/3)+1,3):
			try:
				DIM_5_AA += translate_DNA_codon(DIM_5_seq[i:i+3])
			except:
				print("fail",DIM_5_seq[i:i+3])
		if reference_AA != DIM_5_AA:
			dim5_synonymy = 0
			snp_list = FEATURE[2]
			non_synonymous_snp_set = [snp_list,reference_AA,DIM_5_AA,DIM_5_seq]
			non_syn_feats[2] = non_synonymous_snp_set

		# assign feature to one of four possible groups, based on synonymy
		#	of both mutant sequence and the DIM5 sequence
		if mutant_synonymy == 'T' and dim5_synonymy == 'T':
			non_syn_snps[n][0].append(non_syn_feats)
		elif mutant_synonymy == 0 and dim5_synonymy == 'T':
			non_syn_snps[n][1].append(non_syn_feats)
		# if the DIM_5 strain and mutant strain are both non-synonymous with the reference, yet also non-synonymous with eachother, then
		#	include in the new mutation set
		elif (mutant_synonymy == 0 and dim5_synonymy == 0 and (non_syn_feats[2][1] != non_syn_feats[2][2])):
			non_syn_snps[n][1].append(non_syn_feats)
#			print(non_syn_feats[2][1],non_syn_feats[2][2],"are not equivalent",sep="\n")
		# if the DIM_5 strain and mutant strain are both non-synonymous with the reference, yet still synonymous with eachother
		#	then this non-synonymous CDS is not the cause of the mutation
		elif (mutant_synonymy == 0 and dim5_synonymy == 0 and (non_syn_feats[2][1] == non_syn_feats[2][2])):
			non_syn_snps[n][2].append(non_syn_feats)
#			print(non_syn_feats[2][1],non_syn_feats[2][2],"are equivalent",sep="\n")
		else:
			non_syn_snps[n][3].append(non_syn_feats)
			


#########################################################################################
#	results and output
#########################################################################################
synonymous = 0
partial = 0
non = 0
revert = 0
for contig in non_syn_snps:
	synonymous += len(contig[0])	
	partial += len(contig[1])	
	non += len(contig[2])	
	revert += len(contig[3])
print()
print("RESULTS")
print("synonymous with reference :",synonymous,"new mutation :",partial,"same mutation as DIM_5 parent :",non,"revert from DIM_5 mutation back to reference :",revert, sep="\t")
print()		

with open(OUT_FILE,'w') as out:
	for n in range(0,7):
		out.write("Supercontig 12."+str(n+1)+"\n")
		if (len(non_syn_snps[n][1]) > 0):
			out.write("NEW_MUTATION"+"\n")
			new_mutation_set = non_syn_snps[n][1]
			for feature in new_mutation_set:
				insertion_count = 0
				deletion_count = 0
				snp_count = 0
				feature_name = feature[0][0]
				snp_list = feature[1][0]
				ref_AA = feature[1][1]
				mutant_AA = feature[1][2]
				out.write(feature_name+"\n")
				for snp in snp_list:
					out.write("\t".join(map(str,snp)) + "\n")
					if (len(snp[1]) == len(snp[2])):
						snp_count += 1
					elif (len(snp[1]) < len(snp[2])):
						insertion_count += 1
					if (len(snp[1]) > len(snp[2])):
						deletion_count += 1
				for line in (difflib.Differ().compare(ref_AA.splitlines(1),mutant_AA.splitlines(1))):
					out.write(line.strip('\n')+"\n")
				if ("-" in mutant_AA):
					out.write('snps:'+"\t"+str(snp_count)+"\t"+'insertions:'+"\t"+str(insertion_count)+"\t"+'deletions:'+"\t"+str(deletion_count)+"\t"+'STOP'+"\t"+"\n")
				else:
					out.write('snps:'+"\t"+str(snp_count)+"\t"+'insertions:'+"\t"+str(insertion_count)+"\t"+'deletions:'+"\t"+str(deletion_count)+"\t"+"\n")
				out.write("\n")
		else:
			out.write('null'+"\n")

out.close()

print()
print('mpileup reference and user calculated reference sequence agreement:')
print('matches',matches)
print('conflicts',conflicts)
print(str(conflicts/(matches+conflicts)*100),'% error rate')

sys.exit()
