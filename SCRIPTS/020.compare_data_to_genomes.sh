#!/bin/bash

# Run all scripts in order
cd `pwd`
BASE="$( dirname "$( dirname "$( echo `pwd` )" )" )"
FILES="$BASE/ESSENTIAL/SCRIPTS"

# run scripts as
#	$FILES/{ script_name_here }

# either
# $FILES/021.single_end_map_to_OR.sh
# or
$FILES/021.paired_end_map_to_OR.sh
$FILES/022.bed_filtered_bam.sh
$FILES/024.find_snps_between_genomes.sh
$FILES/025.parse_snps.sh
