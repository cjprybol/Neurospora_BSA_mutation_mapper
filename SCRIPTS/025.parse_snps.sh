#!/bin/bash

cd `pwd`
BASE="$( dirname "$( dirname "$( echo `pwd` )" )" )"
FILES="$BASE"/VCF_OUTPUT/*.vcf


DIR="$BASE"/VCF_OUTPUT

# find snps that are present in both parents (i.e. are divergent from the reference genome in the same way)
#	and only output 'CHROM', 'POS', 'ID', 'REF', 'ALT'


CONFIG_FILE="$BASE/ESSENTIAL/config.txt"

oak_ridge="$( grep "^OR:" $CONFIG_FILE | perl -pe 's/OR://' )"
mauriceville="$( grep "^MV:" $CONFIG_FILE | perl -pe 's/MV://' )"

cat "$DIR/$oak_ridge.bed_filtered.vcf" "$DIR/$mauriceville.bed_filtered.vcf" | gawk '{FS="\t"}{OFS="\t"}{if ($1 ~ /^Supercontig_12.[1-7]$/) print $1,$2,$3,$4,$5}' | sort | uniq -d > "$DIR/both_parents.snps"


for f in "$DIR/$mauriceville.bed_filtered.vcf"
do
	in_file=${f##*/}

	cleaned_name=$( echo "$f" | perl -pe 's/\.vcf$/\.cleaned\.vcf/' )

	# remove ## info header lines from vcf file
	sed -r -i '/^##/d' $f
	
	# remove all snps not in supercontigs 1 thru 7
	gawk '{FS="\t"}{OFS="\t"}{if ($1 ~ /^Supercontig_12.[1-7]$/) print $0}' $f > $f.contig_filtered

	# remove snps present in both parents that are divergent from the reference genome
	# python script works by locating two lines in a row that are repeated (done by inserting the
	#	the snps in order) and removing both
	cat $f.contig_filtered "$DIR/both_parents.snps" | sort -k 1,3 > $f.combined

	# python3 025.1.clean_snps.py [file with duplicate snps inserted] [out file with snps common in both parents removed]
	python3 025.1.clean_snps.py "$f.combined" "$cleaned_name"

	# clean up vcf file for easier use in later analysis
	python3 025.2.parse_snps.py "$cleaned_name"

	echo "##fileformat=VCFv4.2" > "$DIR/temp"
	head -1 $f >> "$DIR/temp"
	cat "$cleaned_name" >> "$DIR/temp"
	mv "$DIR/temp" "$cleaned_name"

	# clean up temp files
	rm "$DIR/both_parents.snps"
	rm $f.contig_filtered
	rm $f.combined

done
