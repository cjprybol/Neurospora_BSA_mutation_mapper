#bin/bash

# Author: Cameron Prybol
# Created: 2014.06.25
# Last Updated:
# Description: 

#BASE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd | perl -pe 's|/[A-Z_]+/[A-Z_]+$||g')"
BASE="/lustre1/escratch1/cprybol1_Jul_30"
FILES="$BASE"/OR_MAP_BAM/*

#########################################################
#	if output folder doesn't exist, make it
########################################################

if [ ! -d "$BASE/GFF_FILTERED_BAM_ORIGINAL" ];
        then
                mkdir "$BASE/GFF_FILTERED_BAM_ORIGINAL"
                echo "> created directory $BASE/GFF_FILTERED_BAM_ORIGINAL"
fi

##############################################################
#	assign directory variable names
##############################################################

IN_DIR="$BASE/OR_MAP_BAM"

# determine number of files, used for creating readable output to user
num_files=$(ls -1 "$IN_DIR" | wc -l)

i=1

OUT_DIR="$BASE/GFF_FILTERED_BAM_ORIGINAL"

###############################################################################################################
#	
##############################################################################################################

for f in $FILES
do
        echo "> Processing file $i of $num_files" 

	#remove the path prefix on the file name
	in_file=${f##*/}
	out_file=$(echo "$in_file" | sed -e 's/\.sorted.*/\.gff_filtered/')

	#create a temporary .sam file, with un-mapped reads removed
	samtools view -F 4 "$IN_DIR/$in_file" > "$OUT_DIR/$i.tmp.sam"

	#run python script to remove all reads not in gff features
	python3 05.sub_script.py "$BASE/ESSENTIAL/REF_GENOMES/Ncrassa_OakRidge/neurospora_crassa_or74a_12_transcripts.gff3" \
	"$OUT_DIR/$i.tmp.sam" "$OUT_DIR/$i.filtered.sam"
	
	rm "$OUT_DIR/$i.tmp.sam"

	#convert filtered .sam to sorted .bam
	samtools view -bT "$BASE/ESSENTIAL/REF_GENOMES/Ncrassa_OakRidge/neurospora_crassa_or74a_12_supercontigs.fasta" \
	"$OUT_DIR/$i.filtered.sam" | samtools sort - "$OUT_DIR/$out_file"

	rm "$OUT_DIR/$i.filtered.sam"

        let i=i+1

done
