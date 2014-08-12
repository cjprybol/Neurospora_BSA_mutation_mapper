#bin/bash

# Author: Cameron Prybol
# Created: 2014.06.17
# Last Updated: 2014.06.19
# Description: maps fastQ sequence files of N. crassa to the oak ridge genome

#BASE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd | perl -pe 's|/[A-Z_]+/[A-Z_]+$||g')"
BASE="/lustre1/escratch1/cprybol1_Jul_30"
FILES="$BASE"/ESSENTIAL/MERGED_FASTQ/*

##############################################################
#	build index from reference genome
############################################################

/usr/local/bwa/latest/bwa index "$BASE/ESSENTIAL/REF_GENOMES/Ncrassa_OakRidge/neurospora_crassa_or74a_12_supercontigs.fasta"


###############################################################
#	if bam output folder doesn't exist, make it
###############################################################

if [ ! -d "$BASE/OR_MAP_BAM" ];
        then
                mkdir "$BASE/OR_MAP_BAM"
                echo "> created directory $BASE/OR_MAP_BAM"
fi

##############################################################
#	map fasta files to bowtie reference genome files
##############################################################

# determine number of files, used for creating readable output to user
num_files=$(ls -1 "$BASE"/ESSENTIAL/MERGED_FASTQ/ | wc -l)

i=1

BAM_DIR="$BASE/OR_MAP_BAM"

for f in $FILES
do

	# /usr/local/bwa/latest/bwa mem [options] <idxbase> <in1.fq> [in2.fq]
	file=${f##*/}
	file_name=$(echo "$file" | sed -e 's/\.f\(q\|astq\)//')
	echo "> Processing file $i of $num_files : $file"
	/usr/local/bwa/latest/bwa mem "$BASE/ESSENTIAL/REF_GENOMES/Ncrassa_OakRidge/neurospora_crassa_or74a_12_supercontigs.fasta" "$f" | "$BASE/ESSENTIAL/xa2multi.pl" | samtools view -buS - | samtools sort - "$BAM_DIR/$file_name.sorted"
	let i=i+1
done
