#!/bin/bash

# 

BASE="/escratch4/cprybol1/cprybol1_Nov_19"
FILES="$BASE"/VCF_OUTPUT/*

for f in $FILES
do

        # create variable with filename without full directory path
        in_file=${f##*/}

	# create variable with output filename where input bam file extension is swapped for output vcf file extension
        out_vcf=$(echo "$in_file" | sed -e 's/\.bam/\.vcf/')

	### figure out how to run python script from here

done
