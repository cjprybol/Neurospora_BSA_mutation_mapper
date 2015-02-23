#!/usr/local/bin/python3

import sys
import subprocess
import pandas as pd

TRANSCRIPTS = sys.argv[1]
FASTA = sys.argv[2]
OUTFILE = sys.argv[3]

transcripts = pd.io.parsers.read_table(TRANSCRIPTS, index_col=0, header=None)
transcripts['sequence'] = ''

with open(FASTA) as f:
	match = 0
	current_transcript = ""
	sequence = ""
	for line in f:
		line = line.strip()
		if (line.startswith(">")):
			
			# if length > 0, then there must be a full sequence, add to dataframe
			if (len(sequence) > 0):
				transcripts['sequence'][current_transcript] = sequence

			line = line.split("|")[0].strip()[1:]
			if (line in transcripts.index):
				match = 1
				current_transcript = line
				sequence = ""
			else:
				match = 0
				current_transcript = ""
				sequence = ""
		elif (match == 1):
			sequence += line

transcripts.to_csv(OUTFILE, sep="\t")
