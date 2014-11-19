#/bin/bash

# Author: Cameron Prybol
# Created: 2014.07.25
# Last Updated:
# Description: create .bcf and .vcf files from .bam file

#BASE="/lustre1/escratch1/cprybol1_Jul_23/DIM5_SUPPRESSOR_PROJECT"
#BASE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd | perl -pe 's|/[A-Z_]+/[A-Z_]+$||g')"
BASE="/escratch3/cprybol1/cprybol1_Sep_11"
FILES="$BASE"/VCF_OUTPUT/lane3*.filtered

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
num_files=$(ls -1 "$BASE/VCF_OUTPUT" | grep "^lane3*" |  grep "filtered" | wc -l)

i=1

IN_DIR="$BASE/VCF_OUTPUT"
SAM_DIR="$BASE/KB_REGION_FILTERED_BAM"
OUT_DIR="$BASE/SYNONYMY"
REF_GENOME_DIR="$BASE/ESSENTIAL/REF_GENOMES/Ncrassa_OakRidge"
GFF="$REF_GENOME_DIR/neurospora_crassa_or74a_12_transcripts.gff3"
TRANSCRIPT_FILE="$REF_GENOME_DIR/neurospora_crassa_or74a_12_supercontigs.fasta"

for f in $FILES
do
        echo "> Processing file $i of $num_files"   

        #remove the path prefix on the file name
        in_file=${f##*/}
	file_head=$(echo "$in_file" | sed -e 's/\.vcf\.filtered//')
	out_file_base="$file_head.synonymy"


	echo "matching VCF locations to genes and translating for synonymy..."
	python3 08.synonymy.py "$GFF" "$TRANSCRIPT_FILE" "$IN_DIR/$file_head.unique" "$OUT_DIR/$out_file_base"

	# parse out SA chimeric reads
	samtools view "$SAM_DIR/$file_head.kb_region_filtered.sorted.bam" | grep 'SA:Z:' > "$OUT_DIR/$file_head.tmp.sam"

	filter_list="$BASE/ESSENTIAL/FILTER_SITES/$file_head.filter_sites"

	echo "extracting SA chimeric reads in windows of interest"
	python3 08.SA_chimeric_read_mapping.py "$GFF" "$OUT_DIR/$file_head.tmp.sam" "$filter_list" "$OUT_DIR/$out_file_base"

	rm "$OUT_DIR/$file_head.tmp.sam"


        let i=i+1
done

