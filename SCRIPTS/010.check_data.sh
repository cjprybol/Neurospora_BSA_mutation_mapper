#!/bin/bash

# Run all scripts in order
cd `pwd`
BASE="$( dirname "$( dirname "$( echo `pwd` )" )" )"
FILES="$BASE/ESSENTIAL/SCRIPTS"

# run scripts as
#	$FILES/{ script_name_here }

sh $FILES/011.fastqc_reads_pre.sh
sh $FILES/012.lighter_paired_end.sh
#sh $FILES/012.lighter_single_end.sh
sh $FILES/013.fastqc_reads_post.sh
