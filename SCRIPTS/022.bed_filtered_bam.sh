#!/bin/bash

cd `pwd`
BASE="$( dirname "$( dirname "$( echo `pwd` )" )" )"
FILES="$BASE"/OR_MAP_BAM/*.bam

IN_DIR="$BASE/OR_MAP_BAM"
OUT_DIR="$BASE/BED_FILTERED_BAM"
#########################################################
#	if output folder doesn't exist, make it
########################################################

if [ ! -d "$OUT_DIR" ];
        then
                mkdir "$OUT_DIR"
                echo "> created directory $OUT_DIR"
fi

##############################################################
#	assign directory variable names
##############################################################

CONFIG_FILE="$BASE/ESSENTIAL/config.txt"
read_type="$( grep "^read-type:" $CONFIG_FILE | perl -pe 's/read-type://' )"


for f in $FILES
do
	#remove the path prefix on the file name
	in_file=${f##*/}
	out_file=$(echo "$in_file" | sed -e 's/\.sorted.*/\.bed_filtered/')

	# -v	Only report those entries in A that have _no overlaps_ with B.  - Similar to "grep -v" (an homage).
	# -abam	The A input file is in BAM format.  Output will be BAM as well.
	# -b <bed/gff/vcf>
	# i.e. find all reads (in -abam) that are not in -b, and write those to file (filter out regions in -b)
	bedtools intersect -v -abam "$f" -b "$BASE/ESSENTIAL/HETEROCHROMATIN/Sorted_S1_H2K9me3_rseg_default_peaks.bed" > "$OUT_DIR/$out_file.bam.tmp"
	# filter out unmapped reads

	if [ "$read_type" == "se" ]; then

		samtools view -bh -F 1796 -q 20 "$OUT_DIR/$out_file.bam.tmp" | samtools sort - "$OUT_DIR/$out_file"

	elif [ "$read_type" == "pe" ]; then

		samtools view -bh -F 1804 -q 20 "$OUT_DIR/$out_file.bam.tmp" | samtools sort - "$OUT_DIR/$out_file"

	else
		echo "read-type incorrectly specified in $BASE/ESSENTIAL/config.txt"
	fi

	rm "$OUT_DIR/$out_file.bam.tmp"

done
