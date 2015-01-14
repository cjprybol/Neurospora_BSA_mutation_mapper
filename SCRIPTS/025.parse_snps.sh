#!/bin/bash

BASE="/escratch4/cprybol1/cprybol1_Nov_19"
FILES="$BASE"/VCF_OUTPUT/*.gz

for f in $FILES
do

        # create variable with filename without full directory path
        in_file=${f##*/}

	### figure out how to run python script from here
	python3 025.parse_snps.py $f "$BASE"/VCF_OUTPUT/parsed_snps.out

done
