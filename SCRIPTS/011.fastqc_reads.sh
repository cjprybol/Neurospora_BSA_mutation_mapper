#!/bin/bash

# run FASTQC on reads to check quality of data
cd `pwd`
BASE="/escratch4/cprybol1/cprybol1_Jan_21"
FILES="$BASE"/ESSENTIAL/FASTQ/*.fastq


# if output folder does not exist, create it
if [ ! -d "$BASE/FASTQC_OUT" ];
        then
                mkdir "$BASE/FASTQC_OUT"
                echo "> created directory $BASE/FASTQC_OUT"
fi


OUT_DIR="$BASE/FASTQC_OUT"

# run fastqc on all fastq files
for f in $FILES
do

	echo -e "running FASTQC on \n\t $f"
        file=${f##*/}
        /usr/local/fastqc/latest/fastqc --outdir $OUT_DIR $f
	
done
