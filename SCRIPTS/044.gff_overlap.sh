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

	# if snps overlap features in GFF file, report them
	/usr/local/bedtools/latest/bin/bedtools intersect -wo -a $GFF -b $f > "$OUT_DIR/$file_head.all"

	# if anything falls in a CDS window, then it must be checked for amino acid synonymy, so report which transcripts need to be checked
	grep 'CDS' "$OUT_DIR/$file_head.all" | awk '{print $9}' | perl -pe 's/.*?;Parent=(NCU[0-9]+T[0-9])/$1/' > "$OUT_DIR/$file_head.transcript_list"

	# make another file listing all snp info associated with the transcript
	echo -e "CDS\tpos\tref\talt" > "$OUT_DIR/$file_head.snp_list"
	grep 'CDS' "$OUT_DIR/$file_head.all" | awk '{OFS="\t"}{print $9, $11, $13, $14}' | perl -pe 's/ID=CDS:(NCU\d+T\d:\d);\w+=\w+/$1/' >> "$OUT_DIR/$file_head.snp_list"

	# create blank file, overwriting file if it already exists
	> "$OUT_DIR/$file_head.cds_windows"
	while read p; do
		grep "$p" "$GFF" | grep "CDS" | awk '{OFS="\t"}{print $1, $3, $4, $5, $7, $8, $9}' >> "$OUT_DIR/$file_head.cds_windows"
	done < "$OUT_DIR/$file_head.transcript_list"

	/usr/local/bedtools/latest/bin/bedtools intersect -v -a $f -b $GFF > "$OUT_DIR/$file_head.snps_not_in_genes"

done
