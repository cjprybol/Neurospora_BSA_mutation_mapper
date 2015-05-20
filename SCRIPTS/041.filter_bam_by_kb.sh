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


CONFIG_FILE="$BASE/ESSENTIAL/config.txt"

oak_ridge="$( grep "^OR:" $CONFIG_FILE | perl -pe 's/OR://' )"
mauriceville="$( grep "^MV:" $CONFIG_FILE | perl -pe 's/MV://' )"

declare -a FILES=($(ls -1 $IN_DIR/*.bam))


# drop oak ridge and mauricveille files from list of files
index=0
for keyword in ${FILES[@]}; do
        if [ "$keyword" == "$IN_DIR/$oak_ridge.bed_filtered.bam" ] || [ "$keyword" == "$IN_DIR/$mauriceville.bed_filtered.bam" ]; then
                unset FILES[$index]
        fi
        let index++
done


for f in ${FILES[@]}
do

	in_file=${f##*/}
	file_head=$(echo "$in_file" | sed -e 's/\.bed_filtered\.bam//')

	# assign filter list to variable
	filter_list="$BASE/ESSENTIAL/FILTER_SITES/$file_head.filter_sites.bed"

	oak_ridge_file="$IN_DIR/$oak_ridge.bed_filtered.bam"
	mauriceville_file="$IN_DIR/$mauriceville.bed_filtered.bam"

	samtools view -hbu -L $filter_list $f | samtools sort - "$OUT_DIR/$file_head.kb_filtered"
	samtools view -hbu -L $filter_list $oak_ridge_file | samtools sort - "$OUT_DIR/$file_head.OR"
	samtools view -hbu -L $filter_list $mauriceville_file | samtools sort - "$OUT_DIR/$file_head.MV"
	
done
