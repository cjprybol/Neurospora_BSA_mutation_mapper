#/bin/bash

# Author: Cameron Prybol
# Created: 2014.07.25
# Last Updated:
# Description: create .bcf and .vcf files from .bam file

#BASE="/lustre1/escratch1/cprybol1_Jul_23/DIM5_SUPPRESSOR_PROJECT"
#BASE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd | perl -pe 's|/[A-Z_]+/[A-Z_]+$||g')"
BASE="/lustre1/escratch1/cprybol1_Aug_22"
FILES="$BASE"/VCF_OUTPUT/lane3*

#########################################################
#       
########################################################

if [ ! -d "$BASE/SYNONYMY" ];
        then
                mkdir "$BASE/SYNONYMY"
                echo "> created directory $BASE/SYNONYMY"
fi

################################################################
#
################################################################

# determine number of files, used for creating readable output to user
num_files=$(ls -1 "$BASE/VCF_OUTPUT" | grep "^lane3" |  wc -l)

i=1

IN_DIR="$BASE/VCF_OUTPUT"
OUT_DIR="$BASE/SYNONYMY"
REF_GENOME_DIR="$BASE/ESSENTIAL/REF_GENOMES/Ncrassa_OakRidge"
GFF="$REF_GENOME_DIR/neurospora_crassa_or74a_12_transcripts.gff3"
TRANSCRIPT_FILE="$REF_GENOME_DIR/neurospora_crassa_or74a_12_supercontigs.fasta"

for f in $FILES
do
        echo "> Processing file $i of $num_files"   

        #remove the path prefix on the file name
        in_file=${f##*/}
	DIM_5_file="R1.barcode10.$in_file"
	file_head=$(echo "$in_file" | sed -e 's/\.vcf//')
	out_file="$file_head.synonymy.out"

	python3 08.synonymy.py "$GFF" "$TRANSCRIPT_FILE" "$IN_DIR/$in_file" "$IN_DIR/$DIM_5_file" "$OUT_DIR/$out_file"


        let i=i+1
done

