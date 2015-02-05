#!/bin/bash

BASE="/escratch4/cprybol1/cprybol1_Jan_21"
FILES="$BASE"/SNP_MAPPING/*snp_map\.out

OUT_DIR="$BASE/SNP_MAPPING/PARSED_SNP_INFO"

if [ ! -d "$OUT_DIR" ];
        then
                mkdir "$OUT_DIR"
                echo "> created directory $OUT_DIR"
fi

for f in $FILES
do

	# create variable containing filename but without full directory path
	file=${f##*/}

	out_base=$(echo "$file" | perl -pe 's/\.bed_filtered\.snp_map\.out//')

	window=5
	slide=1

	out_file="$out_base.$window""kb_window"

	echo -e "CONTIG\tPOS\tREF\tALT\tMIS" > "$OUT_DIR/$out_base.temp"
	grep 'CONTIG' $f | awk '{OFS="\t"}{print $2,$4,$6,$8,$10}' | sort -k1n,1 -k2n,2 >> "$OUT_DIR/$out_base.temp" 

	python3 033.bucket_count.py "$OUT_DIR/$out_base.temp" "$OUT_DIR/$out_file" $window $slide

	rm "$OUT_DIR/$out_base.temp"

done
