#!/bin/bash

# Run all scripts in order
cd `pwd`
BASE="$( dirname "$( dirname "$( echo `pwd` )" )" )"
FILES="$BASE/ESSENTIAL/SCRIPTS"

CONFIG_FILE="$BASE/ESSENTIAL/config.txt"
read_type="$( grep "^read-type:" $CONFIG_FILE | perl -pe 's/read-type://' )"

# run scripts as
#	$FILES/{ script_name_here }

if [ "$read_type" == "se" ]; then
	echo "single end"
	$FILES/011.fastqc_reads_pre.sh &
	$FILES/012.lighter_single_end.sh &
	wait
elif [ "$read_type" == "pe" ]; then
	echo "paired end"
	$FILES/011.fastqc_reads_pre.sh &
	$FILES/012.lighter_paired_end.sh &
	wait
else
	echo "read-type incorrectly specified in $BASE/ESSENTIAL/config.txt"
fi

$FILES/013.fastqc_reads_post.sh
