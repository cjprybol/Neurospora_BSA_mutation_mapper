#!/bin/bash

# Run all scripts in order
cd `pwd`
BASE="$( dirname "$( dirname "$( echo `pwd` )" )" )"
FILES="$BASE/ESSENTIAL/SCRIPTS"

# run scripts as
#	$FILES/{ script_name_here }

$FILES/041.filter_bam_by_kb.sh
$FILES/042.find_snps_between_genomes.sh
$FILES/043.parse_snps.sh
$FILES/044.gff_overlap.sh
$FILES/045.translate_cds.sh
