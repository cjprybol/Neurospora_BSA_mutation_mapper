#/bin/bash

# Author: Cameron Prybol
# Created: 2014.07.25
# Last Updated:
# Description: create .bcf and .vcf files from .bam file

#BASE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd | perl -pe 's|/[A-Z_]+/[A-Z_]+$||g')"
BASE="/escratch3/cprybol1/cprybol1_Sep_11"
FILES="$BASE"/VCF_OUTPUT/lane3*.filtered

i=1

VCF_DIR="$BASE/VCF_OUTPUT"

for f in $FILES
do
	in_file=${f##*/}
	out_file=$(echo "$in_file" | sed -r 's/\.vcf\.filtered/\.unique/g')
	dim_5_file="R1.barcode10.$in_file"
	
	grep "^[^#]" "$VCF_DIR/$in_file" | awk '{OFS="\t"} {print $1, $2, $4, $5}' > "$VCF_DIR/$i.1.tmp"
	grep "^[^#]" "$VCF_DIR/$dim_5_file" | awk '{OFS="\t"} {print $1, $2, $4, $5}' > "$VCF_DIR/$i.2.tmp"
	gawk 'NR==FNR{a[$0];next}!($0 in a)' "$VCF_DIR/$i.2.tmp" "$VCF_DIR/$i.1.tmp" > "$VCF_DIR/$i.uniq"
	rm "$VCF_DIR/$i.1.tmp"
	rm "$VCF_DIR/$i.2.tmp"
	grep "^[^#]" "$VCF_DIR/$in_file" > "$VCF_DIR/$i.header_less"
	cat "$VCF_DIR/$i.uniq" "$VCF_DIR/$i.header_less" | sort -k1,2 > "$VCF_DIR/$i.sorted_with_duplicates"
	rm "$VCF_DIR/$i.uniq"
	rm "$VCF_DIR/$i.header_less"
	python3 "$BASE/ESSENTIAL/SCRIPTS/07.1.subscript.py" "$VCF_DIR/$i.sorted_with_duplicates" "$VCF_DIR/$out_file"
	rm "$VCF_DIR/$i.sorted_with_duplicates"
	

        let i=i+1

done
