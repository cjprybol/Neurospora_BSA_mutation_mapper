#!/bin/bash

# Run all scripts in order
cd `pwd`
BASE="$( dirname "$( dirname "$( echo `pwd` )" )" )"
FILES="$BASE/ESSENTIAL/SCRIPTS"

echo $FILES

# run scripts as
#	$FILES/{ script_name_here }

$FILES/010.check_data.sh
$FILES/020.compare_data_to_genomes.sh
$FILES/030.determine_snp_profiles_for_samples.sh
