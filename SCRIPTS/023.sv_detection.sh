#!/bin/bash

cd `pwd`
BASE="$( dirname "$( dirname "$( echo `pwd` )" )" )"

IN_DIR="$BASE/OR_MAP_BAM"
OUT_DIR="$BASE/SV_DETECTION"
#########################################################
#	if output folder doesn't exist, make it
########################################################

if [ ! -d "$OUT_DIR" ];
        then
                mkdir "$OUT_DIR"
                echo "> created directory $OUT_DIR"
fi

CONFIG_FILE="$BASE/ESSENTIAL/config.txt"

oak_ridge="$( grep "^OR:" $CONFIG_FILE | perl -pe 's/OR://' )"
mauriceville="$( grep "^MV:" $CONFIG_FILE | perl -pe 's/MV://' )"

declare -a FILES=($(ls -1 $IN_DIR/*.bam))

# drop oak ridge and mauricveille files from list of files
index=0
for keyword in ${FILES[@]}; do
        if [ "$keyword" == "$IN_DIR/$oak_ridge.sorted.bam" ] || [ "$keyword" == "$IN_DIR/$mauriceville.sorted.bam" ]; then
                unset FILES[$index]
        fi
        let index++
done

BAM_FILES=$IN_DIR/*.bam
for f in $BAM_FILES
do

	samtools index $f

done


for f in ${FILES[@]}
do

	in_file=${f##*/}
	oak_ridge_bam="$IN_DIR/$oak_ridge.sorted.bam"

	bam2cfg.pl -g -h tumor.bam normal.bam
	bam2cfg.pl "$oak_ridge_bam" "$f" > "$OUT_DIR/$in_file.cfg"
	breakdancer_max -g  "$OUT_DIR/$in_file.bed" "$OUT_DIR/$in_file.cfg" > "$OUT_DIR/$in_file.breakdancer_out"

	# create pindel config file
	perl -e 'print reverse <>' "$OUT_DIR/$in_file.cfg" > "$OUT_DIR/$in_file.pindel_config"
	awk '{OFS="\t"}{print $3,$9}' "$OUT_DIR/$in_file.pindel_config" | perl -pe 's/^map://g' | perl -pe 's/mean://g' > "$OUT_DIR/tmp"

	input_sample_name=$( echo $in_file | perl -pe 's/_pool.*//g' )
	
	i=1
	while read p; do
		if [ "$i" == 1 ]; then
			echo -e "$p\t$input_sample_name" > "$OUT_DIR/$in_file.pindel_config"
		else
			echo -e "$p\tDim_5" >> "$OUT_DIR/$in_file.pindel_config"
		fi
		i=$((i + 1))
	done <"$OUT_DIR/tmp"

	rm "$OUT_DIR/tmp"

#	pindel -f <reference.fa> -i bam_configuration_file -c <chromosome_name> -o <prefix_for_output_file>
	pindel -f "$BASE/ESSENTIAL/REF_GENOMES/Ncrassa_OakRidge/neurospora_crassa_or74a_12_supercontigs.fasta" -i "$OUT_DIR/$in_file.pindel_config" -c ALL -o "$OUT_DIR/$input_sample_name"

done
