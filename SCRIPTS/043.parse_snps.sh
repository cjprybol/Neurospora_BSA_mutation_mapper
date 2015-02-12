#!/bin/bash

cd `pwd`
BASE="$(dirname "$( dirname "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" )" )"
FILES="$BASE"/VCF_OUTPUT/*.vcf

DIR="$BASE"/VCF_OUTPUT

# find snps that are present in both parents (i.e. are divergent from the reference genome in the same way)
#	and only output 'CHROM', 'POS', 'ID', 'REF', 'ALT'

cat "$DIR/1_dim-5_49-19.bed_filtered.vcf" "$DIR/2_Mauriceville.bed_filtered.vcf" | gawk '{FS="\t"}{OFS="\t"}{if ($1 ~ /^Supercontig_12.[1-7]$/) print $1,$2,$3,$4,$5}' | sort | uniq -d > "$DIR/both_parents.snps"

for f in "$DIR/2_Mauriceville.bed_filtered.vcf"
do
	in_file=${f##*/}

	# remove ## info header lines from vcf file
	sed -r -i '/^##/d' $f
	
	# remove all snps not in supercontigs 1 thru 7
	gawk '{FS="\t"}{OFS="\t"}{if ($1 ~ /^Supercontig_12.[1-7]$/) print $0}' $f > $f.contig_filtered

	# remove snps present in both parents that are divergent from the reference genome
	# python script works by locating two lines in a row that are repeated (done by inserting the
	#	the snps in order) and removing both
	cat $f.contig_filtered "$DIR/both_parents.snps" | sort -k 1,3 > $f.combined

	# python3 025.1.clean_snps.py [file with duplicate snps inserted] [out file with snps common in both parents removed]
	python3 025.1.clean_snps.py $f.combined $f.cleaned
	
	# clean up vcf file for easier use in later analysis
	python3 025.2.parse_snps.py $f.cleaned

	# clean up temp files
	rm $f.contig_filtered
	rm $f.combined
	rm $f.cleaned

done
