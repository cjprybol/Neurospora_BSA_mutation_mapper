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

# determine number of files, used for creating readable output to user
num_files=$(ls -1 "$IN_DIR" | wc -l)

i=1

OUT_DIR="$BASE/BED_FILTERED_BAM"

###############################################################################################################
#	
##############################################################################################################

#for f in $FILES
for f in "$IN_DIR"/4_CPX22_pool_method_1.sorted.bam "$IN_DIR"/5_CPX22_pool_method_2.sorted.bam
do
        echo "> Processing file $i of $num_files" 

	#remove the path prefix on the file name
	in_file=${f##*/}
	out_file=$(echo "$in_file" | sed -e 's/\.sorted.*/\.bed_filtered/')

	/usr/local/bedtools/latest/bin/bedtools intersect -v -abam "$f" -b "$BASE/ESSENTIAL/HETEROCHROMATIN/Sorted_S1_H2K9me3_rseg_default_peaks.bed" > "$OUT_DIR/$out_file.bam.tmp"
	# filter out unmapped reads
	samtools view -bh -F 4 -q 10 "$OUT_DIR/$out_file.bam.tmp" | samtools sort - "$OUT_DIR/$out_file"

	rm "$OUT_DIR/$out_file.bam.tmp"

        let i=i+1

done
