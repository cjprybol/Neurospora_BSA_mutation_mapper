#!/bin/bash

cd `pwd`
BASE="$(dirname "$( dirname "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" )" )"

IN_DIR="$BASE"/GFF_OVERLAP

FILES="$IN_DIR"/*.CDS


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
	file_head=$( echo "$in_file" | perl -pe 's/\.CDS//')

	perl -pe 's/.*?Parent=(NCU[0-9]+T[0-9]).*/$1/' $f > "$OUT_DIR/$file_head.transcripts.temp"

	python3 045.1.extract_fasta_seq.py "$OUT_DIR/$file_head.transcripts.temp" "$TRANSCRIPTS" "$OUT_DIR/$file_head.fasta_extract"
	python3 045.2.translate.py "$f" "$OUT_DIR/$file_head.fasta_extract" "$OUT_DIR/$file_head.translated_CDS"

	rm "$OUT_DIR/$file_head.transcripts.temp"
	rm "$OUT_DIR/$file_head.fasta_extract"


done
