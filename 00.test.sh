#bin/bash

# Author: Cameron Prybol
# Created: 2014.06.17
# Last Updated: 
# Description: 

#BASE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd | perl -pe 's|/[A-Z_]+/[A-Z_]+$||g')"
BASE="/lustre1/escratch1/cprybol1_Jul_30"
FILES="$BASE/ESSENTIAL/SCRIPTS"

#"$FILES/01.map_to_OR.sh" 
"$FILES/02.map_to_MV.sh"
"$FILES/03.remove_MV_snps.sh"
"$FILES/04.mpileup.sh"
"$FILES/05.gff_filtered_bam.sh"
#"$FILES/06.filter_sam_by_kb.sh"
#"$FILES/07.bcf_vcf"
#"$FILES/08.synonymy.sh"
