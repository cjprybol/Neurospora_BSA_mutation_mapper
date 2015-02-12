#!/bin/bash

# map fastq reads to the oak ridge genome
cd `pwd`
BASE="$(dirname "$( dirname "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" )" )"

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

	/usr/local/bowtie2/latest/bin/bowtie2-build "$BASE/ESSENTIAL/REF_GENOMES/Ncrassa_OakRidge/neurospora_crassa_or74a_12_supercontigs.fasta" "$INDEX_DIR/or_index"

##############################################################
#	map fasta files to bwa reference genome files
##############################################################


# only want to grab the R1 files
FILES="$BASE"/ESSENTIAL/FASTQ/*R1.fastq

for f in $FILES
do

	# create variable containing filename but without full directory path
	file=${f##*/}
	folder=$(dirname $f)

	# get file base to append R2.fastq to, allowing for variable handling of the R1/R2 forward/reverse filenames
	base=$(echo "$file" | perl -pe 's/(.*)?\.R1.*/$1/')

	# get base filename without the arbitrary numerical identifier
	clean_base=$(echo "$file" | perl -pe 's/(.*)?-.*/$1/')

	# get reverse file
	rev="$base.R2.fastq"

	# get output file
	out="$clean_base"

	echo "in1: $f"
	echo "in2: $folder/$rev"

	# set max distance for paired end to 3000 !! set this to data !! based on quality control data from sequencing center
	/usr/local/bowtie2/latest/bin/bowtie2 -X 3000 -x "$INDEX_DIR/or_index" -1 $f -2 "$folder/$rev" | samtools view -buS - | samtools sort - "$BAM_DIR/$out.sorted"

done
