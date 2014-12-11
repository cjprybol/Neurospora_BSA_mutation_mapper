#bin/bash

# Author: Cameron Prybol
# Created: 2014.06.17
# Last Updated: 
# Description: 

#BASE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd | perl -pe 's|/[A-Z_]+/[A-Z_]+$||g')"
BASE="/escratch4/cprybol1/cprybol1_Nov_19"
FILES="$BASE/ESSENTIAL/SCRIPTS"

#"$FILES/01.map_to_OR.sh"
#"$FILES/02.1.simulate_MV_reads.sh"
"$FILES/02.2.map_sim_reads_to_OR.sh"
"$FILES/02.3.vcf.sh"
"$FILES/02.4.parse_snps.py"

#"$FILES/03.remove_MV_snps.sh"
#"$FILES/04.bed_filtered_bam.sh"
#"$FILES/05.bedfiltered_mpileup.sh"
#"$FILES/06.filter_sam_by_kb.sh"
#"$FILES/07.bcf_vcf"
#"$FILES/07.1.uniq_snps.sh"
#"$FILES/08.synonymy.sh"
