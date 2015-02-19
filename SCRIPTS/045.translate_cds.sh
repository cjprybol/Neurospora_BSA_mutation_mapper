#!/bin/bash

cd `pwd`
BASE="$( dirname "$( dirname "$( echo `pwd` )" )" )"

IN_DIR="$BASE"/GFF_OVERLAP

FILES="$IN_DIR"/*.cds_windows


#########################################################
#       if output folder doesn't exist, make it
########################################################
OUT_DIR="$IN_DIR/TRANSLATE_CDS"

if [ ! -d "$OUT_DIR" ];
        then
                mkdir "$OUT_DIR"
                echo "> created directory $OUT_DIR"
fi

TRANSCRIPTS="$BASE/ESSENTIAL/REF_GENOMES/Ncrassa_OakRidge/neurospora_crassa_or74a_12_transcripts.fasta"

for f in $FILES
do
	in_file="${f##*/}"
	file_head=$( echo "$in_file" | perl -pe 's/\.cds_windows//')

	cds_file="$f"
	transcript_list="$file_head.transcript_list"
	snps="$file_head.snp_list"

	python3 045.1.extract_fasta_seq.py "$IN_DIR/$transcript_list" "$TRANSCRIPTS" "$OUT_DIR/$file_head.fasta_extract"
	python3 045.2.translate.py "$cds_file" "$IN_DIR/$snps" "$OUT_DIR/$file_head.fasta_extract" "$OUT_DIR/$file_head.translated_CDS"

	rm "$OUT_DIR/$file_head.fasta_extract"


done
