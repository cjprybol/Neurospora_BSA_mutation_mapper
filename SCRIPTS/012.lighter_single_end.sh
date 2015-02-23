#!/bin/bash

# run FASTQC on reads to check quality of data
cd `pwd`
BASE="$( dirname "$( dirname "$( echo `pwd` )" )" )"

FILES="$(ls "$BASE"/ESSENTIAL/FASTQ/*.fastq)"

OUT_DIR="$BASE/LIGHTER_FASTQ"

# if output folder does not exist, create it
if [ ! -d "$OUT_DIR" ];
        then
                mkdir "$OUT_DIR"
                echo "> created directory $OUT_DIR"
fi

# write blank file to append output data to
> "$OUT_DIR/lighter.out"


# run ligther on all fastq files
for f in $FILES
do

	
	file=${f##*/}
	echo -e "running Lighter on \n\t $file" >> "$OUT_DIR/lighter.out"

	# run lighter using the forward and reverse read pairs, with a k-mer size
	# of 19 (see Figure 5 of paper)
	# -K kmer_length genom_size: in this case, the genome size should be relative accurate
	# -od set out directory
	# -trim: allow trimming at ends of low-quality
	# -stable: sequentialize the sampling stage, output the same result with different runs (default: false)
	# -t num_of_threads: number of threads to use (default: 1)
	lighter -r "$f" -K 19 41020000 -od "$OUT_DIR" -trim -stable -t 4 2>> "$OUT_DIR/lighter.out"
	
done
