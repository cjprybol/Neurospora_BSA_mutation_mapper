#!/bin/bash

# map fastq reads to the oak ridge genome

cd `pwd`
BASE="$(dirname "$( dirname "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" )" )"

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

for f in $IN_DIR/3_CPX3_pool.bed_filtered.bam $IN_DIR/4_CPX22_pool_method_1.bed_filtered.bam $IN_DIR/5_CPX22_pool_method_2.bed_filtered.bam
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
	python3 032.map_snps_for_samples.py "$OUT_DIR/$sam_temp" "$OUT_DIR/$out_file" "$BASE/VCF_OUTPUT/2_Mauriceville.parsed_snps.out"

#	rm "$BASE/OR_MAP_BAM_CLEANED/$sam_temp"

done
