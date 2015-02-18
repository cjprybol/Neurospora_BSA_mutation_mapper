#!/bin/bash

# Run all scripts in order
cd `pwd`
BASE="$(dirname	"$( dirname "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" )" )"
FILES="$BASE/ESSENTIAL/SCRIPTS"

echo $FILES

# run scripts as
#	$FILES/{ script_name_here }

sh $FILES/010.check_data.sh
sh $FILES/020.compare_data_to_genomes.sh
sh $FILES/030.determine_snp_profiles_for_samples.sh
sh $FILES/040.check_for_synonymy.sh
