#!/bin/bash

cd `pwd`
BASE="$( dirname "$( dirname "$( echo `pwd` )" )" )"
IN_DIR="$BASE/SNP_MAPPING/PARSED_SNP_INFO"
FILES="$IN_DIR/*.snp_counts"

OUT_DIR="$IN_DIR/GRAPHS"

if [ ! -d "$OUT_DIR" ];
        then
                mkdir "$OUT_DIR"
                echo "> created directory $OUT_DIR"
fi

for f in $FILES
do

	# create variable containing filename but without full directory path
	file=${f##*/}

	file_base=$( echo $file | perl -pe 's/\.snp_counts//' )

	R --vanilla --slave "--args $f $OUT_DIR/ $file_base" < "$BASE/ESSENTIAL/SCRIPTS/034.graph_snps.R"

done
