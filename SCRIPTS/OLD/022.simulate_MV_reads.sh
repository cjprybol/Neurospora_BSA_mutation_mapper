#!/bin/bash
# simulate mauriceville reads from the genome

cd `pwd`
BASE="/escratch4/cprybol1/cprybol1_Jan_21"
FILES="$BASE"/ESSENTIAL/MERGED_FASTQ/*

###############################################################
#	if output folder doesn't exist, make it
###############################################################

if [ ! -d "$BASE/MV_SIM_READS" ];
        then
                mkdir "$BASE/MV_SIM_READS"
                echo "> created directory $BASE/MV_SIM_READS"
fi


##############################################################
#	build index from reference genome
############################################################

# -N INT        number of read pairs [1000000]		10x standard
# -e FLOAT      base error rate [0.020]			no errors
# -r FLOAT      rate of mutations [0.0010]		no mutations
# -R FLOAT      fraction of indels [0.15]		no indels
# -X FLOAT      probability an indel is extended [0.30]	no indels
# -s INT        standard deviation [50]			no std dev?

/usr/local/samtools/latest/wgsim -N10000000 -e0 -r0 -R0 -X0 -s0 "$BASE/ESSENTIAL/REF_GENOMES/Ncrassa_Mauriceville/neurospora_crassa_MV_contigs.fasta" "$BASE/MV_SIM_READS/mv_sim.fq" /dev/null
