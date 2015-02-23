#!/bin/bash

# Run all scripts in order
cd `pwd`
BASE="$( dirname "$( dirname "$( echo `pwd` )" )" )"
FILES="$BASE/ESSENTIAL/SCRIPTS"

echo $FILES

# run scripts as
#	$FILES/{ script_name_here }

sh $FILES/040.check_for_synonymy.sh
