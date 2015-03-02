#!/bin/bash

cd `pwd`
BASE="$( dirname "$( dirname "$( echo `pwd` )" )" )"

###############################################################
#	if output folder doesn't exist, make it
###############################################################

OUT_DIR="$BASE/SNP_MAPPING"

if [ ! -d "$OUT_DIR" ];
        then
                mkdir "$OUT_DIR"
                echo "> created directory $OUT_DIR"
fi

IN_DIR="$BASE"/BED_FILTERED_BAM

PARENT_FILE="$BASE/ESSENTIAL/FASTQ/parents.txt"

oak_ridge="$( grep "^OR:" $PARENT_FILE | perl -pe 's/OR://' )"
mauriceville="$( grep "^MV:" $PARENT_FILE | perl -pe 's/MV://' )"

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

	# create variable containing filename but without full directory path
	file=${f##*/}

	# create variable with output filename
	file_base=$(echo "$file" | sed -e 's/\.bam//')
	sam_temp="$file_base.sam"
	out_file="$file_base.snp_map.out"

	# create temp sam file for python script with pre-cleaned fields for efficiency
	#	view without header |
	#	if reads are on supercontig_12.1 thru 7, return the read_id, supercontig, position, and sequence |
	#	drop 'Supercontig_12.' and leave only the trailing number >
	#	save to file
	echo "writing sam"
	samtools view $f | awk '{OFS="\t"}{ if ($3 ~ /^Supercontig_12.[1-7]$/) print $1,$3,$4,$10}' | sed 's/Supercontig_12\.//'> "$OUT_DIR/$sam_temp"

	echo "cleaning"
	# python3 032.map_snps_for_samples.py [sam_file] [out_file] [snp_list]
	python3 032.map_snps_for_samples.py "$OUT_DIR/$sam_temp" "$OUT_DIR/$out_file" "$BASE/VCF_OUTPUT/$mauriceville.parsed_snps.out"

	rm "$OUT_DIR/$sam_temp"

done
