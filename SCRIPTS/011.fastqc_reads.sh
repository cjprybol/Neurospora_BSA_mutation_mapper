#!/bin/bash

# run FASTQC on reads to check quality of data
cd `pwd`
BASE="$(dirname "$( dirname "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" )" )"

FILES="$(ls "$BASE"/ESSENTIAL/FASTQ/*.fastq)"

OUT_DIR="$BASE/FASTQC_OUT"

# if output folder does not exist, create it
if [ ! -d "$OUT_DIR" ];
        then
                mkdir "$OUT_DIR"
                echo "> created directory $OUT_DIR"
fi



# run fastqc on all fastq files
for f in $FILES
do

	echo -e "running FASTQC on \n\t $f"
	file=${f##*/}
	/usr/local/fastqc/latest/fastqc --outdir $OUT_DIR $f
	
done
