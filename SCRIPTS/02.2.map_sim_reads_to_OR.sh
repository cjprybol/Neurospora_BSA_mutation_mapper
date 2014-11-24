#bin/bash

# Author: Cameron Prybol
# Created: 2014.06.17
# Last Updated: 2014.06.19
# Description: maps fastQ sequence files of N. crassa to the oak ridge genome

#BASE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd | perl -pe 's|/[A-Z_]+/[A-Z_]+$||g')"
BASE="/escratch4/cprybol1/cprybol1_Nov_19"
FILES="$BASE"/MV_SIM_READS/mv_sim.fq

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

BAM_DIR="$BASE/OR_MAP_BAM"

for f in $FILES
do

	# /usr/local/bwa/latest/bwa mem [options] <idxbase> <in1.fq> [in2.fq]
	file=${f##*/}
	file_name=$(echo "$file" | sed -e 's/\.f\(q\|astq\)//')
	echo "> Processing file $i of $num_files : $file"
	/usr/local/bwa/latest/bwa mem -a "$BASE/ESSENTIAL/REF_GENOMES/Ncrassa_OakRidge/neurospora_crassa_or74a_12_supercontigs.fasta" "$f" | samtools view -buS - | samtools sort - "$BAM_DIR/$file_name.sorted"
	let i=i+1
done
