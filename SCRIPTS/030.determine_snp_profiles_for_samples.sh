#!/bin/bash

# Run all scripts in order
cd `pwd`
BASE="$( dirname "$( dirname "$( echo `pwd` )" )" )"
FILES="$BASE/ESSENTIAL/SCRIPTS"

# run scripts as
#	$FILES/{ script_name_here }

$FILES/032.map_snps_for_samples.sh
$FILES/033.parse_snp_data.sh
