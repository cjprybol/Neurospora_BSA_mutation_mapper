#!/bin/bash

# Run all scripts in order
cd `pwd`
BASE="$( dirname "$( dirname "$( echo `pwd` )" )" )"
FILES="$BASE/ESSENTIAL/SCRIPTS"

CONFIG_FILE="$BASE/ESSENTIAL/config.txt"
read_type="$( grep "^read-type:" $CONFIG_FILE | perl -pe 's/read-type://' )"

# run scripts as
#       $FILES/{ script_name_here }

if [ "$read_type" == "se" ]; then
        echo "single end"
	$FILES/021.single_end_map_to_OR.sh
elif [ "$read_type" == "pe" ]; then
        echo "paired end"
	$FILES/021.paired_end_map_to_OR.sh
else
        echo "read-type incorrectly specified in $BASE/ESSENTIAL/config.txt"
fi

$FILES/022.bed_filtered_bam.sh
$FILES/024.find_snps_between_genomes.sh
$FILES/025.parse_snps.sh
