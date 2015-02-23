#!/bin/bash

# map fastq reads to the oak ridge genome
cd `pwd`
BASE="$( dirname "$( dirname "$( echo `pwd` )" )" )"

INDEX_DIR="$BASE/OR_INDEX"
BAM_DIR="$BASE/OR_MAP_BAM"

###############################################################
#	if output folder for index doesn't exist, make it
###############################################################

if [ ! -d "$INDEX_DIR" ];
        then
                mkdir "$INDEX_DIR"
                echo "> created directory $INDEX_DIR"
fi

###############################################################
#	if output folder doesn't exist, make it
###############################################################

if [ ! -d "$BAM_DIR" ];
        then
                mkdir "$BAM_DIR"
                echo "> created directory $BAM_DIR"
fi

##############################################################
#	build index from reference genome
############################################################

	bowtie2-build "$BASE/ESSENTIAL/REF_GENOMES/Ncrassa_OakRidge/neurospora_crassa_or74a_12_supercontigs.fasta" "$INDEX_DIR/or_index"

##############################################################
#	map fasta files to bwa reference genome files
##############################################################


FILES="$(ls "$BASE"/LIGHTER_FASTQ/*.fastq)"

for f in $FILES
do


	# create variable containing filename but without full directory path
	file=${f##*/}

	# get base name without extension
	out=$(echo "$file" | perl -pe 's/\.fastq//')

	bowtie2 -p 4 -x "$BASE/OR_INDEX/or_index" -1 $f | samtools view -buS - | samtools sort - "$BAM_DIR/$out.sorted"

done
