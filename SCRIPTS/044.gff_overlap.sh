#!/bin/bash

cd `pwd`
BASE="$(dirname "$( dirname "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" )" )"

FILES="$BASE"/KB_FILTERED_VCF_OUTPUT/CLEANED/*.vcf

IN_DIR="$BASE"/KB_FILTERED_VCF_OUTPUT/CLEANED


#########################################################
#       if output folder doesn't exist, make it
########################################################
OUT_DIR="$BASE/GFF_OVERLAP"

if [ ! -d "$OUT_DIR" ];
        then
                mkdir "$OUT_DIR"
                echo "> created directory $OUT_DIR"
fi

GFF="$BASE/ESSENTIAL/REF_GENOMES/Ncrassa_OakRidge/neurospora_crassa_or74a_12_transcripts.gff3"

for f in $FILES
do
	in_file="${f##*/}"
	file_head="$( echo $in_file | perl -pe 's/\.cleaned\.vcf//')"
	out_all="$file_head.all"
	out_CDS="$file_head.CDS"
	out_no_overlap="$file_head.no_overlap"

	/usr/local/bedtools/latest/bin/bedtools intersect -wo -a $GFF -b $f > "$OUT_DIR/$out_all"
	grep 'CDS' "$OUT_DIR/$out_all" > "$OUT_DIR/$out_CDS"
	/usr/local/bedtools/latest/bin/bedtools intersect -v -a $f -b $GFF > "$OUT_DIR/$out_no_overlap"

done
