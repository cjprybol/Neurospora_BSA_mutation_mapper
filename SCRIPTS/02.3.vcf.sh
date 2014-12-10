#/bin/bash

# Author: Cameron Prybol
# Created: 2014.07.25
# Last Updated:
# Description: create .bcf and .vcf files from .bam file

#BASE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd | perl -pe 's|/[A-Z_]+/[A-Z_]+$||g')"
BASE="/escratch4/cprybol1/cprybol1_Nov_19"
FILES="$BASE"/MV_MAP_BAM/*.bam

#########################################################
#       if vcf output folder doesn't exist, make it
########################################################

if [ ! -d "$BASE/VCF_OUTPUT" ];
        then
                mkdir "$BASE/VCF_OUTPUT"
                echo "> created directory $BASE/VCF_OUTPUT"
fi

################################################################
#
################################################################

# determine number of files, used for creating readable output to user
num_files=$(ls -1 "$BASE/MV_MAP_BAM" | grep "bam$" |  wc -l)

i=1

IN_DIR="$BASE/MV_MAP_BAM"
OUT_VCF_DIR="$BASE/VCF_OUTPUT"
REF_GENOME_DIR="$BASE/ESSENTIAL/REF_GENOMES/Ncrassa_OakRidge"

for f in $FILES
do
        echo "> Processing file $i of $num_files"   

        #remove the path prefix on the file name
        in_file=${f##*/}
        out_vcf=$(echo "$in_file" | sed -e 's/\.bam/\.vcf/')

#		#	mpilup options
#		#	-u	Generate uncompressed VCF/BCF output, which is preferred for piping
#		#	-g	Compute genotype likelihoods and output them in the binary call format (BCF). As of v1.0, this is BCF2 which is incompatible with the BCF1 format produced by previous (0.1.x) versions of samtools.
#		#	-f	reference fasta	
#		#	bcftools call options
#		#	-v	output variant sites only
#		#	-m	multiallelic caller, alternative model for multiallelic and rare−variant calling designed to overcome known limitations in −c calling model
#		#	−O v	output vcf format
#		#	-o	output filename
#	
#		samtools mpileup -ugf "$REF_GENOME_DIR/neurospora_crassa_or74a_12_supercontigs.fasta" "$IN_DIR/$in_file" | bcftools call -vmO z -o "$OUT_VCF_DIR/$out_vcf.gz"
#		tabix -p vcf "$OUT_VCF_DIR/$out_vcf.gz"
#		bcftools stats -F "$REF_GENOME_DIR/neurospora_crassa_or74a_12_supercontigs.fasta" -s - "$OUT_VCF_DIR/$out_vcf.gz" > "$OUT_VCF_DIR/$out_vcf.gz.stats"
#	
	gunzip < "$OUT_VCF_DIR/$out_vcf.gz" | grep -v '^##' > "$OUT_VCF_DIR/$out_vcf.temp"

        let i=i+1
done
