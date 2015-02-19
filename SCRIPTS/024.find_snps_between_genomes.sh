#!/bin/bash

cd `pwd`
BASE="$( dirname "$( dirname "$( echo `pwd` )" )" )"

OUT_VCF_DIR="$BASE/VCF_OUTPUT"

#########################################################
#       if output folder doesn't exist, make it
########################################################

if [ ! -d "$OUT_VCF_DIR" ];
        then
                mkdir "$OUT_VCF_DIR"
                echo "> created directory $OUT_VCF_DIR"
fi

REF_GENOME_DIR="$BASE/ESSENTIAL/REF_GENOMES/Ncrassa_OakRidge"
IN_DIR="$BASE/BED_FILTERED_BAM"
FILES="$BASE"/BED_FILTERED_BAM/*.bam

for f in $FILES
do

        # create variable with filename without full directory path
        in_file=${f##*/}

	# create variable with output filename where input bam file extension is swapped for output vcf file extension
        out_vcf=$(echo "$in_file" | sed -e 's/\.bam/\.vcf/')

#		#	mpilup options
#		#	-u	Generate uncompressed VCF/BCF output, which is preferred for piping
#		#	-g	Compute genotype likelihoods and output them in the binary call format (BCF). As of v1.0, this is BCF2 which is incompatible with the BCF1 format produced by previous (0.1.x) versions of samtools.
#		#	-f	reference fasta	
#		#	bcftools call options
#		#	−O v	output vcf format
#		#	-o	output filename
#	
	samtools mpileup -ugf "$REF_GENOME_DIR/neurospora_crassa_or74a_12_supercontigs.fasta" "$IN_DIR/$in_file" | bcftools call -vmO z -o "$OUT_VCF_DIR/$out_vcf.gz"
	tabix -p vcf "$OUT_VCF_DIR/$out_vcf.gz"
	bcftools stats -F "$REF_GENOME_DIR/neurospora_crassa_or74a_12_supercontigs.fasta" -s - "$OUT_VCF_DIR/$out_vcf.gz" > "$OUT_VCF_DIR/$out_vcf.gz.stats"
	
	gunzip "$OUT_VCF_DIR/$out_vcf.gz"
	rm "$OUT_VCF_DIR/$out_vcf.gz.tbi"

done
