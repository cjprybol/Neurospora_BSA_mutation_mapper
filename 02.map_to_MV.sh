#bin/bash

# Author: Cameron Prybol
# Created: 2014.06.17
# Last Updated: 2014.06.19
# Description: maps fastQ sequence files of N. crassa to the Mauriceville genome

#BASE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd | perl -pe 's|/[A-Z_]+/[A-Z_]+$||g')"
BASE="/lustre1/escratch1/cprybol1_Jul_30"
FILES="$BASE"/ESSENTIAL/MERGED_FASTQ/*

#########################################################
#	if mauriceville genome index folder doesn't exist, make it
########################################################

if [ ! -d "$BASE/MV_BOWTIE_INDEX" ];
	then
		mkdir "$BASE/MV_BOWTIE_INDEX"
		echo "> created directory $BASE/MV_BOWTIE_INDEX"
fi

##############################################################
#	build bowtie index from reference genome
############################################################

bowtie2-build "$BASE/ESSENTIAL/REF_GENOMES/Ncrassa_Mauriceville/neurospora_crassa_MV_contigs.fasta" "$BASE/MV_BOWTIE_INDEX/mv_index"


###############################################################
#	if bam output folder doesn't exist, make it
###############################################################

if [ ! -d "$BASE/MV_MAP_BAM" ];
        then
                mkdir "$BASE/MV_MAP_BAM"
                echo "> created directory $BASE/MV_MAP_BAM"
fi


##############################################################
#	map fasta files to bowtie reference genome files
##############################################################


# determine number of files, used for creating readable output to user
num_files=$(ls -1 "$BASE"/ESSENTIAL/MERGED_FASTQ/ | wc -l)

i=1

INDEX_DIR="$BASE/MV_BOWTIE_INDEX"

BAM_DIR="$BASE/MV_MAP_BAM"

for f in $FILES
do
	file=${f##*/}
	file_name=$(echo "$file" | sed -e 's/\.f\(q\|astq\)//')
	echo "> Processing file $i of $num_files : $file"
	let i=i+1
#	using standard settings because we aren't trying to map any mutation sequences,
#	only 100% matches i.e. regular sensitivity and we run in end-to-end mode
	bowtie2 -x "$INDEX_DIR/mv_index" -U $f | samtools view -bS - | samtools sort - "$BAM_DIR/$file_name.sorted"
done
