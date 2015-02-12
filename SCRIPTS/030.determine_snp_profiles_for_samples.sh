#!/bin/bash

# Run all scripts in order
cd `pwd`
BASE="$(dirname "$( dirname "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" )" )"
FILES="$BASE/ESSENTIAL/SCRIPTS"

# run scripts as
#	$FILES/{ script_name_here }

sh $FILES/031.clean_bam_files.sh
