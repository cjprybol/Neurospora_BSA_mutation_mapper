#!/bin/bash

cd `pwd`
BASE="$( dirname "$( dirname "$( echo `pwd` )" )" )"
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

	echo -e "CONTIG\tPOS\tREF\tALT\tMIS" > "$OUT_DIR/$out_base.snp_counts"

	grep '^CONTIG' $f | awk '{FS="\t"}{OFS="\t"}{print $2,$4,$6,$8,$10}' >> "$OUT_DIR/$out_base.snp_counts"

done
