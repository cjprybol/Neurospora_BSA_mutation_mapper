#!/bin/bash

# map simulated mauriceville reads to the oak ridge genome

BASE="/escratch4/cprybol1/cprybol1_Nov_19"
FILES="$BASE"/MV_SIM_READS/mv_sim.fq

##############################################################
#	build index from reference genome
############################################################

	/usr/local/bwa/latest/bwa index "$BASE/ESSENTIAL/REF_GENOMES/Ncrassa_OakRidge/neurospora_crassa_or74a_12_supercontigs.fasta"


###############################################################
#	if output folder doesn't exist, make it
###############################################################

if [ ! -d "$BASE/MV_SIM_MAP_BAM" ];
        then
                mkdir "$BASE/MV_SIM_MAP_BAM"
                echo "> created directory $BASE/MV_SIM_MAP_BAM"
fi

##############################################################
#	map fasta files to bowtie reference genome files
##############################################################

BAM_DIR="$BASE/MV_SIM_MAP_BAM"

for f in $FILES
do

	# create variable containing filename but without full directory path
	file=${f##*/}

	# create variable with same name but without .fastq and/or .fq file extension
	file_name=$(echo "$file" | sed -e 's/\.f\(q\|astq\)//')

	# -a            output all alignments for SE or unpaired PE
	/usr/local/bwa/latest/bwa mem "$BASE/ESSENTIAL/REF_GENOMES/Ncrassa_OakRidge/neurospora_crassa_or74a_12_supercontigs.fasta" "$f" | samtools view -buS - | samtools sort - "$BAM_DIR/$file_name.sorted"

done
