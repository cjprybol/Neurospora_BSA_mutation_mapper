#bin/bash

# Author: Cameron Prybol
# Created: 2014.07.17
# Last Updated: 
# Description: creates oak ridge .bam files with mauriceville snps removed

#BASE="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd | perl -pe 's|/[A-Z_]+/[A-Z_]+$||g')"
BASE="/lustre1/escratch1/cprybol1_Jul_30"
FILES="$BASE"/MV_MAP_BAM/*

#########################################################
#	if output folder doesn't exist, make it
########################################################

if [ ! -d "$BASE/OR_MINUS_MV_SNPS_BAM" ];
        then
                mkdir "$BASE/OR_MINUS_MV_SNPS_BAM"
                echo "> created directory $BASE/OR_MINUS_MV_SNPS_BAM"
fi


MV_DIR="$BASE/MV_MAP_BAM"
OR_DIR="$BASE/OR_MAP_BAM"

# determine number of files, used for creating readable output to user
num_files=$(ls -1 "$MV_DIR" | wc -l)

i=1

OUT_DIR="$BASE/OR_MINUS_MV_SNPS_BAM"

for f in $FILES
do
        echo "> Processing fileset $i of $num_files" 

	#remove the path prefix on the file name
	in_file=${f##*/}


##############################################################################################################
#	view mauriceville .bam file ('-F 4' remove unmapped reads '-q 10' remove anything mapped with less than
#	90% confidence of location) | get lines of reads that are exact matches | extract read ID |
#	sort | output to temp file
##############################################################################################################
	echo "> extracting read identifiers from exact_hits to Mauriceville"
	samtools view -F 4 -q 10 "$MV_DIR/$in_file" | grep NM:i:0 | awk '{print $1}' | sort > "$OUT_DIR/$i.tmp.MV.junk"


##############################################################################################################
#	extract read ID's from oak ridge .bam file
##############################################################################################################
        echo "> extracting read identifiers from oak ridge sam file"
        samtools view "$OR_DIR/$in_file" | awk '{print $1}' | sort > "$OUT_DIR/$i.tmp.OR.junk"


##############################################################################################################
#	output contents of both temp files together | sort them | find duplicates (via -d option on uniq) >
#	write to temp file (note this fails if sort is not used before uniq -d)
##############################################################################################################
        echo "> extracting read identifiers from the mauriceville exact hits file also in the Oak_Ridge .bam file"
        cat "$OUT_DIR/$i.tmp.OR.junk" "$OUT_DIR/$i.tmp.MV.junk" | sort | uniq -d > "$OUT_DIR/$i.duplicates"

        # remove first set of temp files
        rm "$OUT_DIR/$i.tmp.OR.junk"
        rm "$OUT_DIR/$i.tmp.MV.junk"

	# create sam file needed for next step
	samtools view "$OR_DIR/$in_file" > "$OUT_DIR/$i.tmp.sam"


##############################################################################################################
#	print out full oak ridge .sam file and duplicates
#	sort (-T to specify out directory, meant to solve size issue encountered with temp files made by sort command ....
#	... -k1,1 to sort only on column one of each line) > output to temp file
##############################################################################################################
	cat "$OUT_DIR/$i.tmp.sam" "$OUT_DIR/$i.duplicates" | sort -T "$OUT_DIR" -k1,1 > "$OUT_DIR/$i.sorted_with_duplicates.sam"

        # remove second set of temp files
        rm "$OUT_DIR/$i.duplicates"
	rm "$OUT_DIR/$i.tmp.sam"

##############################################################################################################
#	this python script creates a new .sam file without the snps. please see the script
#	for more detail
##############################################################################################################
        echo "> eliminating duplicates into Oak Ridge .sam file"
        python3 03.sub_script.py "$OUT_DIR/$i.sorted_with_duplicates.sam" "$OUT_DIR/$i.MV_hits_removed.sam"
        rm "$OUT_DIR/$i.sorted_with_duplicates.sam"

        file_head=$(echo "$in_file" | sed -e 's/\.sorted\.bam//')
	
##############################################################################################################
#	add header lines back to the .sam file, convert to .bam, and sort
##############################################################################################################
	samtools view -bT "$BASE/ESSENTIAL/REF_GENOMES/Ncrassa_OakRidge/neurospora_crassa_or74a_12_supercontigs.fasta" "$OUT_DIR/$i.MV_hits_removed.sam" | samtools sort - "$OUT_DIR/$file_head.minus_mv_snps.sorted"
	rm "$OUT_DIR/$i.MV_hits_removed.sam"
        let i=i+1

done
