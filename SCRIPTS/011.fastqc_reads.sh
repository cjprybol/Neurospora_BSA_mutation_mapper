#bin/bash

# run FASTQC on reads to check quality of data

BASE="/escratch4/cprybol1/cprybol1_Nov_19"
FILES="$BASE"/ESSENTIAL/MERGED_FASTQ/*

# create output folder to direct results to

if [ ! -d "$BASE/FASTQC_OUT" ];
        then
                mkdir "$BASE/FASTQC_OUT"
                echo "> created directory $BASE/FASTQC_OUT"
fi


OUT_DIR="$BASE/FASTQC_OUT"

for f in $FILES
do

	echo -e "running FASTQC on \n\t $f"
        file=${f##*/}
        /usr/local/fastqc/latest/fastqc --outdir $OUT_DIR $f
	
done
