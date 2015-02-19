#!/bin/bash

cd `pwd`
BASE="$( dirname "$( dirname "$( echo `pwd` )" )" )"

IN_DIR="$BASE/BED_FILTERED_BAM"
OUT_DIR="$BASE/KB_FILTERED_BAM"

#########################################################
#	if output folder doesn't exist, make it
########################################################

if [ ! -d "$OUT_DIR" ];
        then
                mkdir "$OUT_DIR"
                echo "> created directory $OUT_DIR"
fi

### here
for f in "$IN_DIR/3_CPX3_pool.bed_filtered.bam" "$IN_DIR/4_CPX22_pool_method_1.bed_filtered.bam" "$IN_DIR/5_CPX22_pool_method_2.bed_filtered.bam"
do

	in_file=${f##*/}
	file_head=$(echo "$in_file" | sed -e 's/\.bed_filtered\.bam//')

	echo $file_head

	# assign filter list to variable
	filter_list="$BASE/ESSENTIAL/FILTER_SITES/$file_head.filter_sites.bed"

### here
	oak_ridge_file="$IN_DIR/1_dim-5_49-19.bed_filtered.bam"
	mauriceville_file="$IN_DIR/2_Mauriceville.bed_filtered.bam"

	samtools view -hbu -L $filter_list $f | samtools sort - "$OUT_DIR/$file_head.kb_filtered"
	samtools view -hbu -L $filter_list $oak_ridge_file | samtools sort - "$OUT_DIR/$file_head.OR"
	samtools view -hbu -L $filter_list $mauriceville_file | samtools sort - "$OUT_DIR/$file_head.MV"
	
done
