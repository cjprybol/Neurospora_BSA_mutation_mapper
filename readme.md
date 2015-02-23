##DIM5 Analysis
---

### About
This project was created in order to find candidates for causitive mutations that suppress hypersensitivity to DNA damaging agents in the DIM5 knockout strain of Neurospora Crassa

### setup project
clone the github repository into a directory called `ESSENTIAL`

```
mkdir {desired working directory}
cd {desired working directory}
git clone https://github.com/cprybol/DIM5_suppressor_mapping.git ESSENTIAL

```

### required software (all must be in user path)
- Perl 5
	- used `perl 5, version 14, subversion 1 (v5.14.1) built for x86_64-linux-thread-multi`
- Python2
	- must be the default system python (required by Bowtie for mapping reads to genome)
- Python3
	- recommended installation method is to simply get [Anaconda](http://continuum.io/downloads#34)
- pandas for python3
	- included in Anaconda
- BioPython
	- not included in Anaconda
	- `conda install biopython`
- [Samtools](http://www.htslib.org/)
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [Bedtools](http://bedtools.readthedocs.org/en/latest/)
- [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [Lighter](https://github.com/mourisl/Lighter)



### Prepare environment
1. add fastq files to `ESSENTIAL/FASTQ`
	- follow naming conventions specified in `ESSENTIAL/FASTQ/naming_conventions.txt`
2. specify reference parent and divergent parent filenames in the `ESSENTIAL/FASTQ/parents.txt`
3. edit `ESSENTIAL/SCRIPTS/010.check_data.sh` to run proper fastq error correction script
	- this git repo contains scripts for paired end and single end datasets
	- if using paired end data, remove the `#` pound sign comment before `sh $FILES/012.lighter_paired_end.sh`
	- if using single end data, remove the `#` pound sign comment before `sh $FILES/012.lighter_single_end.sh`
	- ensure the script you are not running is commented out with a `#`
3. edit `ESSENTIAL/SCRIPTS/020.compare_data_to_genomes.sh` to run proper mapping script
	- this git repo contains scripts for mapping paired end and single end reads
	- if running paired end script, remove the `#` pound sign comment before `sh $FILES/021.paired_end_map_to_OR.sh`
		- also, edit the min and max fragment size parameters in the `ESSENTIAL/SCRIPTS/021.paired_end_map_to_OR.sh` file for the bowtie2 command (line 65)
	- if running single end script, remove the `#` pound sign comment before `sh $FILES/021.single_end_map_to_OR.sh`
	- ensure the script you are not running is commented out with a `#`

### steps to run
1. run `ESSENTIAL/SCRIPTS/001.master.sh`
2. graph output data from snp mapping located in `SNP_MAPPING/PARSED_SNP_INFO`
	- make a scatter plot of RATIO ~ KB for each CONTIG, locate regions of interest
3. list your regions of interest in .bed format files in the `ESSENTIAL/FILTER_SITES` directory
	- follow guidelines listed in the readme.txt file in that folder
4. run `ESSENTIAL/SCRIPTS/002.master.sh`
	- obtain output in `GFF_OVERLAP` and `GFF_OVERLAP/TRANSLATE_CDS` folders
	- `GFF_OVERLAP/*.all` lists all GFF features that high quality snps overlap
		- format: {full GFF entry}{full vcf entry for snp}
	- `GFF_OVERLAP/*.snps_not_in_genes` lists all snps that fall in euchromatin regions but do not overlap any GFF features
		- these may hit promotors or other factors outside the gene body, and may be of interest if a mutation is not found in the GFF features
	- `GFF_OVERLAP/TRANSLATE_CDS/*.translated_CDS` outputs the translated AA sequence for all snps falling in coding sequences, and shows if the snp produces a non-synonymous output
