#!/usr/local/bin/python3

import sys
import subprocess
import pandas as pd
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import difflib

CDS = sys.argv[1]
SNP = sys.argv[2]
FASTA_SEQ = sys.argv[3]
OUTFILE = sys.argv[4]


# read in CDS data
# contig	feature_type	start	stop	strand	frame	ID
# Supercontig_12.4        CDS     2674344 2675596 -       2       ID=CDS:NCU06389T0:4;Parent=NCU06389T0
cds_data = pd.io.parsers.read_table(CDS, header=None)


# ID=CDS:(NCU06389T0:4);Parent=NCU06389T0 <- grab CDS name between ()
def parse_cds(x):
	return ':'.join(x.split(";")[0].split("=")[1].split(":")[1:])
cds_data['CDS'] = cds_data[6].apply(lambda x: parse_cds(x))


# ID=CDS:NCU06389T0:4;Parent=(NCU06389T0) <- grab transcript name between ()
def parse_transcript(x):
	return x.split(";")[1].split("=")[1]
cds_data['transcript'] = cds_data[6].apply(lambda x: parse_transcript(x))


# drop contig, feature type, and ID (no longer needed)
cds_data.drop([0,1,6], axis=1, inplace=True)


# rename columns
cds_data.columns = ['start', 'stop', 'strand', 'frame', 'CDS','transcript']


# reorder columns
cds_data = cds_data[['transcript', 'CDS', 'start', 'stop', 'strand', 'frame']]


# read in snp_data
# CDS     pos     ref     alt
# NCU06389T0:4    2674491 T       A
snp_data = pd.io.parsers.read_table(SNP, header=0)


# merge snp info with CDS info on CDS as key
transcript_data = pd.merge(cds_data, snp_data, on='CDS', how='left')


# NCU06389T0:(4) <- only keep CDS number on transcript name
def clean_cds(x):
	return x.split(":")[1]
transcript_data['CDS'] = cds_data['CDS'].apply(lambda x: clean_cds(x))



# reset datatype
transcript_data['CDS'] = transcript_data['CDS'].astype(int)


# determine length of CDS for verification later (this was for testing only, not needed for production)
transcript_data['length'] = cds_data.apply(lambda row: abs(row['start'] - row['stop']) + 1, axis=1)


# import fasta sequence for transcripts
fasta_data = pd.io.parsers.read_table(FASTA_SEQ, index_col=0, header=0)


# rename index to match other dataframe
fasta_data.index.names = ['transcript']


# subset to make a reference for quickly calling which strand the transcript is on
transcript_info = transcript_data[transcript_data['CDS'] == 1][['transcript','strand']]
transcript_info.set_index('transcript', inplace=True)
transcript_list = list(transcript_info.index)
reverse_complement_list = list(transcript_info[transcript_info['strand'] == "-"].index)


# snps are in reference to + strand only, so invert - strand transcripts to simplify later steps
for i in reverse_complement_list:
	fasta_data.at[i,'sequence'] = str(Seq(fasta_data.at[i,'sequence']).reverse_complement())

f = open(OUTFILE, 'w')


# iterate over each transcript to manipulate in snps and check for synonymy
for i in transcript_list:

	# subset data for current transcript
	transcript_set = transcript_data[transcript_data['transcript'] == i]

	# offset data to make everything base 0
	transcript_start = transcript_set['start'].min()
	transcript_set['start'] =  transcript_set['start'].apply(lambda x: x - transcript_start)
	transcript_set['stop'] =  transcript_set['stop'].apply(lambda x: x - transcript_start)
	transcript_set['pos'] =  transcript_set['pos'].apply(lambda x: x - transcript_start)

	max_cds = transcript_set['CDS'].max()

	# reverse numerical order of - strand CDS
	if i in reverse_complement_list:
		transcript_set['CDS'] = transcript_set['CDS'].apply(lambda x: abs(x - max_cds) + 1)

	transcript_set.set_index('CDS', inplace=True)


	# if more than one CDS, there are introns that are spliced out between the CDS windows
	# offset the positions of start, stop, and snp position to remove gaps
	if max_cds > 1:
		for j in range(2,max_cds+1):
			prior_stop = transcript_set.at[j-1,'stop']
			offset = transcript_set.at[j,'start'] - prior_stop
			transcript_set.at[j,'start'] = transcript_set.at[j,'start'] - offset + 1
			transcript_set.at[j,'stop'] = transcript_set.at[j,'stop'] - offset + 1
			transcript_set.at[j,'pos'] = transcript_set.at[j,'pos'] - offset + 1


	# parse out the portion of the transcript for each individual CDS window
	transcript_set['sequence'] = ''
	for j in transcript_set.index:
		start = transcript_set.at[j,'start']
		stop = transcript_set.at[j,'stop']
		transcript_set.at[j,'sequence'] = fasta_data.at[i,'sequence'][start:stop+1]
		#print(len(transcript_set.at[j,'sequence']))
	

	# splice in snps to create modified sequence
	transcript_set['modified_sequence'] = ''
	for j in transcript_set.index:
		if pd.notnull(transcript_set.at[j,'pos']):
			pos = int(transcript_set.at[j,'pos']) - int(transcript_set.at[j,'start'])
			seq = transcript_set.at[j,'sequence']
			current = seq[pos:pos+1]
			ref = transcript_set.at[j,'ref']
			alt = transcript_set.at[j,'alt']
			#print(current, ref, alt, sep="\t")
			transcript_set.at[j,'modified_sequence'] = seq[:pos] + alt + seq[pos+len(ref):]
		else:
			transcript_set.at[j,'modified_sequence'] = transcript_set.at[j,'sequence']


	# build original and snp-modified sequences
	original_seq = ''
	new_seq = ''
	for j in transcript_set.index:
		original_seq += transcript_set.at[j,'sequence']
		new_seq +=  transcript_set.at[j,'modified_sequence']
	original_seq = Seq(original_seq)
	new_seq = Seq(new_seq)
	

	# if on the - strand, flip back before translating
	if i in reverse_complement_list:
		original_seq = original_seq.reverse_complement()
		new_seq = new_seq.reverse_complement()

	# translate
	original_seq = str(original_seq.translate())
	new_seq = str(new_seq.translate())

	# split every n  nucleotides to feed into difflib compare
	n = 50
	ref_list = [(str(i+1) + "." + original_seq[i:i+n]) for i in range(0, len(original_seq), n)]
	new_list = [(str(i+1) + "." + new_seq[i:i+n]) for i in range(0, len(new_seq), n)]

	# use difflib to neatly check for changes to the AA sequence
	diff = difflib.Differ().compare(ref_list, new_list)

	# write output to file
	print(i, "\n".join(diff), "\n", sep="\n", file=f)
