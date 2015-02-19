#!/bin/bash

cd `pwd`
BASE="$( dirname "$( dirname "$( echo `pwd` )" )" )"

FILES="$BASE"/KB_FILTERED_VCF_OUTPUT/*kb_filtered.vcf

DIR="$BASE"/KB_FILTERED_VCF_OUTPUT


#########################################################
#       if output folder doesn't exist, make it
########################################################
OUT_DIR="$DIR/CLEANED"

if [ ! -d "$OUT_DIR" ];
        then
                mkdir "$OUT_DIR"
                echo "> created directory $OUT_DIR"
fi


for f in $FILES
do
	in_file="${f##*/}"
	file_head="$( echo $in_file | perl -pe 's/\.kb_filtered\.vcf//')"

	cat "$DIR/$file_head.OR.vcf" "$DIR/$file_head.MV.vcf" | grep "^[^#]" |gawk '{FS="\t"}{OFS="\t"}{print $1,$2,$3,$4,$5}' | sort | uniq > "$DIR/$file_head.either_parent"

	# remove ## info header lines from vcf file
	grep "^[^#]" $f > "$DIR/$in_file.headerless"
	
	# remove snps present in either parent that are divergent from the reference genome
	# python script works by locating two lines in a row that are repeated (done by inserting the
	#	the snps in order) and removing both
	cat "$DIR/$in_file.headerless" "$DIR/$file_head.either_parent" | sort -k 1,3 > "$DIR/$file_head.combined"

	# python3 025.1.clean_snps.py [file with duplicate snps inserted] [out file with snps common in both parents removed]
	python3 043.1.clean_snps.py "$DIR/$file_head.combined" "$OUT_DIR/$file_head.cleaned"
		
	# clean up vcf file for easier use in later analysis
	python3 043.2.parse_snps.py "$OUT_DIR/$file_head.cleaned" "$OUT_DIR/$file_head.parsed_snps.out"

	#reform VCF file from filtered data
	grep "^#" $f > "$OUT_DIR/$file_head.cleaned.vcf"
	cat "$OUT_DIR/$file_head.parsed_snps.out" >> "$OUT_DIR/$file_head.cleaned.vcf"
	
	# clean up temp files
	rm "$DIR/$file_head.either_parent"
	rm "$DIR/$in_file.headerless"
	rm "$DIR/$file_head.combined"
	rm "$OUT_DIR/$file_head.cleaned"
	rm "$OUT_DIR/$file_head.parsed_snps.out"

done
