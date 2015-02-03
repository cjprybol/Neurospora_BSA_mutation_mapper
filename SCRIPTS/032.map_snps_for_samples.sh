#!/bin/bash

# map fastq reads to the oak ridge genome

BASE="/escratch4/cprybol1/cprybol1_Jan_21"
FILES="$BASE"/OR_MAP_BAM_CLEANED/*.bam

###############################################################
#	if output folder doesn't exist, make it
###############################################################

OUT_DIR="$BASE/SNP_MAPPING"

if [ ! -d "$OUT_DIR" ];
        then
                mkdir "$OUT_DIR"
                echo "> created directory $OUT_DIR"
fi


for f in $FILES
do

	# create variable containing filename but without full directory path
	file=${f##*/}

	# create variable with output filename
	file_base=$(echo "$file" | sed -e 's/\.sorted\.cleaned\.bam//')
	sam_temp="$file_base.sam"
	out_file="$file_base.snp_map.out"

	# create temp sam file for python script with pre-cleaned fields for efficiency
	#	view without header |
	#	if reads are on supercontig_12.1 thru 7, return the read_id, supercontig, position, and sequence |
	#	drop 'Supercontig_12.' and leave only the trailing number >
	#	save to file
#	echo "writing sam"
#	samtools view $f | awk '{OFS="\t"}{ if ($3 ~ /^Supercontig_12.[1-7]$/) print $1,$3,$4,$10}' | sed 's/Supercontig_12\.//'> "$BASE/OR_MAP_BAM_CLEANED/$sam_temp"

#	echo "cleaning"
	# python3 032.map_snps_for_samples.py [sam_file] [out_file] [snp_list]
	python3 032.map_snps_for_samples.py "$BASE/OR_MAP_BAM_CLEANED/$sam_temp" "$OUT_DIR/$out_file" "$BASE/VCF_OUTPUT/parsed_snps.out"

#	rm "$BASE/OR_MAP_BAM_CLEANED/$sam_temp"

done
