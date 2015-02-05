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
	echo $f

#		out_file=$(echo "$file" | perl -pe 's/\.snp_map.out_(12.[1-7])/\.snp_data_$1/')
#	
#		echo "$OUT_DIR/$out_file"	
#	
#		echo -e "position\tref_match_count\talt_match_count\tread_mismatch_count" > $OUT_DIR/$out_file
#		grep '^position' $f | awk '{OFS="\t"}{print $2,$4,$6,$8}' >> $OUT_DIR/$out_file

done
