#!/bin/bash

# Run all scripts in order
cd `pwd`
BASE="$(dirname "$( dirname "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" )" )"
FILES="$BASE/ESSENTIAL/SCRIPTS"

# run scripts as
#	$FILES/{ script_name_here }

sh $FILES/041.filter_sam_by_kb.sh
sh $FILES/042.find_snps_between_genomes.sh
sh $FILES/043.parse_snps.sh
sh $FILES/044.gff_overlap.sh
sh $FILES/045.translate_cds.sh
