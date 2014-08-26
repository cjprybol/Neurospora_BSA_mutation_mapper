#bin/bash

# Author: Cameron Prybol
# Created: 2014.07.24
# Last Updated: 
# Description: creates new .sam files based on user specified kb-windows derived from the snp graphs

#BASE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd | perl -pe 's|/[A-Z_]+/[A-Z_]+$||g')"
BASE="/lustre1/escratch1/cprybol1_Aug_22"
FILES="$BASE"/BED_FILTERED_BAM/lane3*

#########################################################
#	if output folder doesn't exist, make it
########################################################

if [ ! -d "$BASE/KB_REGION_FILTERED_BAM" ];
        then
                mkdir "$BASE/KB_REGION_FILTERED_BAM"
                echo "> created directory $BASE/KB_REGION_FILTERED_BAM"
fi


IN_DIR="$BASE/BED_FILTERED_BAM"

# determine number of files, used for creating readable output to user
num_files=$(ls -1 "$IN_DIR" | grep "lane3" | wc -l)

i=1

OUT_DIR="$BASE/KB_REGION_FILTERED_BAM"

for f in $FILES
do


	echo "> Processing fileset $i of $num_files" 

		#lane3-index01-ATCACG-CPx3.bed_filtered.bam
		in_file=${f##*/}
		file_head=$(echo "$in_file" | sed -e 's/\.bed_filtered\.bam//')
		out_file="$file_head.kb_region_filtered"

		#create temporary .sam file to read from
		samtools view "$IN_DIR/$in_file" > "$OUT_DIR/$i.tmp.sam"

		# assign filter list to variable
		filter_list="$BASE/ESSENTIAL/FILTER_SITES/$file_head.filter_sites"

		# filter the sam file
		python3 06.filter.py "$OUT_DIR/$i.tmp.sam" "$filter_list" "$OUT_DIR/$out_file.sam"

		# convert sam file to .bam file
		samtools view -bT "$BASE/ESSENTIAL/REF_GENOMES/Ncrassa_OakRidge/neurospora_crassa_or74a_12_supercontigs.fasta" "$OUT_DIR/$out_file.sam" | samtools sort - "$OUT_DIR/$out_file.sorted"

		rm "$OUT_DIR/$i.tmp.sam"
		rm "$OUT_DIR/$out_file.sam"
	
        let i=i+1

done

##############################################################################################
#	now create a DIM-5 filter set that matches each of the filter sites for each index
#	need to do this because want to compare DIM-5 sites to index filter sites for each index
##############################################################################################

# grab DIM_5.barcode prefix
DIM_5_file=$(ls -1 "$IN_DIR" | grep "R1")
DIM_5_pre=$(echo "$DIM_5_file" | sed -e 's/\.bed_filtered\.bam//')

echo "> Creating DIM_5 .sam file"

#create temporary .sam file to read from
samtools view "$IN_DIR/$DIM_5_file" > "$OUT_DIR/DIM_5.tmp.sam"


i=1

for f in $FILES
do


	echo "> Creating DIM_5 filtered files based on filter sites for each sample: $i of $num_files" 

		#lane3-index01-ATCACG-CPx3.bed_filtered.bam
		ref_file=${f##*/}
		file_head=$(echo "$ref_file" | sed -e 's/\.bed_filtered\.bam//')

		out_file="$DIM_5_pre.$file_head.kb_region_filtered"

		# assign filter list to variable
		filter_list="$BASE/ESSENTIAL/FILTER_SITES/$file_head.filter_sites"

		# filter the sam file
		python3 06.filter.py "$OUT_DIR/DIM_5.tmp.sam" "$filter_list" "$OUT_DIR/$out_file.sam"

		# convert sam file to .bam file
		samtools view -bT "$BASE/ESSENTIAL/REF_GENOMES/Ncrassa_OakRidge/neurospora_crassa_or74a_12_supercontigs.fasta" "$OUT_DIR/$out_file.sam" | samtools sort - "$OUT_DIR/$out_file.sorted"

		rm "$OUT_DIR/$out_file.sam"
	
        let i=i+1

done
rm "$OUT_DIR/DIM_5.tmp.sam"
