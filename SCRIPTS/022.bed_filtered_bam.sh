#bin/bash

BASE="/escratch4/cprybol1/cprybol1_Jan_21"
FILES="$BASE"/OR_MAP_BAM/*.bam

#########################################################
#	if output folder doesn't exist, make it
########################################################

if [ ! -d "$BASE/BED_FILTERED_BAM" ];
        then
                mkdir "$BASE/BED_FILTERED_BAM"
                echo "> created directory $BASE/BED_FILTERED_BAM"
fi

##############################################################
#	assign directory variable names
##############################################################

IN_DIR="$BASE/OR_MAP_BAM"

OUT_DIR="$BASE/BED_FILTERED_BAM"

for f in $FILES
do
	#remove the path prefix on the file name
	in_file=${f##*/}
	out_file=$(echo "$in_file" | sed -e 's/\.sorted.*/\.bed_filtered/')

	# -v	Only report those entries in A that have _no overlaps_ with B.  - Similar to "grep -v" (an homage).
	# -abam	The A input file is in BAM format.  Output will be BAM as well.
	# -b <bed/gff/vcf>
	# i.e. find all reads (in -abam) that are not in -b, and write those to file (filter out regions in -b)
	/usr/local/bedtools/latest/bin/bedtools intersect -v -abam "$f" -b "$BASE/ESSENTIAL/HETEROCHROMATIN/Sorted_S1_H2K9me3_rseg_default_peaks.bed" > "$OUT_DIR/$out_file.bam.tmp"
	# filter out unmapped reads
	samtools view -bh -F 4 -q 10 "$OUT_DIR/$out_file.bam.tmp" | samtools sort - "$OUT_DIR/$out_file"

	rm "$OUT_DIR/$out_file.bam.tmp"

done
