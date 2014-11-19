#bin/bash

# Author: Cameron Prybol
# Created: 2014.06.25
# Last Updated:
# Description: 

#BASE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd | perl -pe 's|/[A-Z_]+/[A-Z_]+$||g')"
BASE="/escratch3/cprybol1/cprybol1_Sep_11"
FILES="$BASE"/OR_MAP_BAM/*

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

# determine number of files, used for creating readable output to user
num_files=$(ls -1 "$IN_DIR" | wc -l)

i=1

OUT_DIR="$BASE/BED_FILTERED_BAM"

###############################################################################################################
#	
##############################################################################################################

for f in $FILES
do
        echo "> Processing file $i of $num_files" 

	#remove the path prefix on the file name
	in_file=${f##*/}
	out_file=$(echo "$in_file" | sed -e 's/\.sorted.*/\.bed_filtered/')

	/usr/local/bedtools/latest/bin/bedtools intersect -v -abam "$f" -b "$BASE/ESSENTIAL/HETEROCHROMATIN/Sorted_S1_H2K9me3_rseg_default_peaks.bed" > "$OUT_DIR/$out_file.bam.tmp"
	# filter out unmapped reads
	samtools view -h -F 4 -b "$OUT_DIR/$out_file.bam.tmp" > "$OUT_DIR/$out_file.bam"
	rm "$OUT_DIR/$out_file.bam.tmp"

        let i=i+1

done


FILES="$BASE"/OR_MINUS_MV_SNPS_BAM/*

#########################################################
#	if output folder doesn't exist, make it
########################################################

if [ ! -d "$BASE/BED_AND_SNP_FILTERED_BAM" ];
        then
                mkdir "$BASE/BED_AND_SNP_FILTERED_BAM"
                echo "> created directory $BASE/BED_AND_SNP_FILTERED_BAM"
fi

##############################################################
#	assign directory variable names
##############################################################

IN_DIR="$BASE/OR_MINUS_MV_SNPS_BAM"

# determine number of files, used for creating readable output to user
num_files=$(ls -1 "$IN_DIR" | wc -l)

i=1

OUT_DIR="$BASE/BED_AND_SNP_FILTERED_BAM"

###############################################################################################################
#	
##############################################################################################################

for f in $FILES
do
        echo "> Processing file $i of $num_files" 

	#remove the path prefix on the file name
	in_file=${f##*/}
	out_file=$(echo "$in_file" | sed -e 's/\.sorted.*/\.bed_filtered/')

	/usr/local/bedtools/latest/bin/bedtools intersect -v -abam "$f" -b "$BASE/ESSENTIAL/HETEROCHROMATIN/Sorted_S1_H2K9me3_rseg_default_peaks.bed" > "$OUT_DIR/$out_file.bam.tmp"
	# filter out unmapped reads
	samtools view -h -F 4 -b "$OUT_DIR/$out_file.bam.tmp" > "$OUT_DIR/$out_file.bam"
	rm "$OUT_DIR/$out_file.bam.tmp"

        let i=i+1

done
