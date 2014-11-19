#bin/bash

# Author: Cameron Prybol
# Created: 2014.08.12
# Description: maps fastQ sequence files of N. crassa to the mauriceville

#BASE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd | perl -pe 's|/[A-Z_]+/[A-Z_]+$||g')"
BASE="/escratch3/cprybol1/cprybol1_Sep_11"
FILES="$BASE"/ESSENTIAL/MERGED_FASTQ/*

##############################################################
#	build index from reference genome
############################################################

/usr/local/bwa/latest/bwa index "$BASE/ESSENTIAL/REF_GENOMES/Ncrassa_Mauriceville/neurospora_crassa_MV_contigs.fasta"


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

BAM_DIR="$BASE/MV_MAP_BAM"

for f in $FILES
do

	# /usr/local/bwa/latest/bwa mem [options] <idxbase> <in1.fq> [in2.fq]
	file=${f##*/}
	file_name=$(echo "$file" | sed -e 's/\.f\(q\|astq\)//')
	echo "> Processing file $i of $num_files : $file"
	/usr/local/bwa/latest/bwa mem "$BASE/ESSENTIAL/REF_GENOMES/Ncrassa_Mauriceville/neurospora_crassa_MV_contigs.fasta" "$f" | samtools view -buS - | samtools sort - "$BAM_DIR/$file_name.sorted"
	let i=i+1
done
