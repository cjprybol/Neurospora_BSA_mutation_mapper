#bin/bash

# Author: Cameron Prybol
# Created: 2014.07.17
# Last Updated: 
# Description: run mpileup on oak ridge and mauriceville .bam files, and create consensus file
#		as well as the count files used to create plots in R

#BASE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd | perl -pe 's|/[A-Z_]+/[A-Z_]+$||g')"
BASE="/lustre1/escratch1/cprybol1_Jul_30"
FILES="$BASE"/OR_MAP_BAM/*

#########################################################
#	if output folder doesn't exist, make it
########################################################

if [ ! -d "$BASE/MPILEUP" ];
        then
                mkdir "$BASE/MPILEUP"
                echo "> created directory $BASE/MPILEUP"
fi


OR_DIR="$BASE/OR_MAP_BAM"
SNPLESS_DIR="$BASE/OR_MINUS_MV_SNPS_BAM"

# determine number of files, used for creating readable output to user
num_files=$(ls -1 "$OR_DIR" | wc -l)

i=1

OUT_BASE="$BASE/MPILEUP"
REF_FASTA="$BASE/ESSENTIAL/REF_GENOMES/Ncrassa_OakRidge/neurospora_crassa_or74a_12_supercontigs.fasta"
GFF="$BASE/ESSENTIAL/REF_GENOMES/Ncrassa_OakRidge/neurospora_crassa_or74a_12_transcripts.gff3"

for f in $FILES
do


	echo "> Processing fileset $i of $num_files" 


	in_file=${f##*/}
	file_head=$(echo "$in_file" | sed -e 's/\.sorted\.bam//')
	or_file=$(ls -1 "$OR_DIR" | grep "$file_head")
	snpless_file=$(ls -1 "$SNPLESS_DIR" | grep "$file_head")

	# account for variable naming for lane3 samples and provided parent strain fastq filenames
	if [[ "$echo $file_head" == lane3* ]];then
		dir=$(echo "$file_head" | perl -pe 's/lane3-index([0-9]+)-[ACGT]{6}-.*/INDEX\1/')
	else
		dir=$(echo "$file_head" | perl -pe 's/(.*?)\.sorted\.bam/\1/' | tr "[a-z]" "[A-Z]")
	fi
	
	OUT_DIR="$OUT_BASE/$dir"
	echo "$OUT_DIR"


	#########################################################
	#	if output folder doesn't exist, make it
	########################################################

	if [ ! -d "$OUT_DIR" ];
		then
			mkdir "$OUT_DIR"
			echo "> created directory $OUT_DIR"
	fi



	#remove the path prefix on the file name
	samtools mpileup -f "$REF_FASTA" "$OR_DIR/$or_file" > "$OUT_DIR/$file_head.original_OR.mpileup"
	samtools mpileup -f "$REF_FASTA" "$SNPLESS_DIR/$snpless_file" > "$OUT_DIR/$file_head.snpless_OR.mpileup"
	
	COUNTER=1
        while [  $COUNTER -lt 8 ]; do
#		echo "$COUNTER"
		awk '$1 == "Supercontig_12.'$COUNTER'" {print $0}' "$OUT_DIR/$file_head.original_OR.mpileup" > "$OUT_DIR/$file_head.original_OR.Supercontig_12.$COUNTER"
		awk '$1 == "Supercontig_12.'$COUNTER'" {print $0}' "$OUT_DIR/$file_head.snpless_OR.mpileup" > "$OUT_DIR/$file_head.snpless_OR.Supercontig_12.$COUNTER"
		let COUNTER=COUNTER+1
	done

	python3 04.1.mpileup_div.py "$OUT_DIR/$file_head.original_OR" "$OUT_DIR/$file_head.snpless_OR" "$OUT_DIR/$file_head.mpilup"
	rm "$OUT_DIR/$file_head.original_OR.Supercontig_12.*"
	rm "$OUT_DIR/$file_head.snpless_OR.Supercontig_12.*"
	rm "$OUT_DIR/$file_head.original_OR.mpileup"
	rm "$OUT_DIR/$file_head.snpless_OR.mpileup"
	python3 04.2.filter_by_gff.py "$GFF" "$OUT_DIR/$file_head.mpilup.Supercontig_12."
	        COUNTER=1
        while [  $COUNTER -lt 8 ]; do
                echo "$COUNTER"
		python3 04.3.bucket_count.py  "$OUT_DIR/$file_head.mpilup.Supercontig_12.$COUNTER"
                let COUNTER=COUNTER+1
        done

        let i=i+1

done
