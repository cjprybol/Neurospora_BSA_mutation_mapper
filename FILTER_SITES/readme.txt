Supercontig_12.1	475000	825000
Supercontig_12.4	1500000	4700000



above this point, you should list the ranges you wish to keep for
	evaluating snps for possible causitive mutations
regions of the genome outside of these windows will be dropped from
	analysis





use the standard bed format of:
chromosome_name		start_of_window_to_keep		end_of_window_to_keep

seperated by 1 tab







name files by replacing the file extension of the original fastq files
	with .filter_sites.bed)

example filter site file names:

	if starting files are single end fastq:
			lane3-index01-ATCACG-CPx3.fastq
		should have the filter bed file named as
			lane3-index01-ATCACG-CPx3.filter_sites.bed

	if starting files are paired end fastq:
			3_CPX3_pool-21583301.R1.fastq
			3_CPX3_pool-21583301.R2.fastq
		should have a single filter bed file named as
			3_CPX3_pool-21583301.filter_sites.bed
