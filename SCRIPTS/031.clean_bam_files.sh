#bin/bash

# map fastq reads to the oak ridge genome

BASE="/escratch4/cprybol1/cprybol1_Nov_19"
FILES="$BASE"/OR_MAP_BAM/*

###############################################################
#	if output folder doesn't exist, make it
###############################################################

if [ ! -d "$BASE/OR_MAP_BAM_CLEANED" ];
        then
                mkdir "$BASE/OR_MAP_BAM_CLEANED"
                echo "> created directory $BASE/OR_MAP_BAM_CLEANED"
fi

BAM_DIR="$BASE/OR_MAP_BAM_CLEANED"

for f in $FILES
do

	# create variable containing filename but without full directory path
	file=${f##*/}

	# create variable with output filename
	file_name=$(echo "$file" | sed -e 's/\.bam/\.cleaned/')

	# view mauriceville .bam file ('-F 4' remove unmapped reads '-q 10' remove anything 
	#	mapped with less than 90% confidence of location)
	samtools view -bh -F 4 -q 10 $f | samtools sort - $BAM_DIR/$file_name

done
