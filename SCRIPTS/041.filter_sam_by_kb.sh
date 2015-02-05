#bin/bash

BASE="/escratch4/cprybol1/cprybol1_Jan_21"

#########################################################
#	if output folder doesn't exist, make it
########################################################

IN_DIR="$BASE/BED_FILTERED_BAM"
OUT_DIR="$BASE/KB_FILTERED_BAM"

if [ ! -d "$OUT_DIR" ];
        then
                mkdir "$OUT_DIR"
                echo "> created directory $OUT_DIR"
fi

#for f in "$IN_DIR/3_CPX3_pool.bed_filtered.bam" "$IN_DIR/4_CPX22_pool_method_1.bed_filtered.bam" "$IN_DIR/5_CPX22_pool_method_2.bed_filtered.bam"
for f in "$IN_DIR/3_CPX3_pool.bed_filtered.bam"
do

	in_file=${f##*/}
	file_head=$(echo "$in_file" | sed -e 's/\.bed_filtered\.bam//')

	# assign filter list to variable
	filter_list="$BASE/ESSENTIAL/FILTER_SITES/$file_head.filter_sites.bed"

	samtools view -h -L $filter_list $f > "$OUT_DIR/$file_head.tmp.sam"

done
