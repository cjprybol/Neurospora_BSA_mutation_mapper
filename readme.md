##DIM5 Analysis
---

This study was supported in part by resources and technical expertise from the Georgia Advanced Computing Resource Center, a partnership between the University of Georgia’s Office of the Vice President for Research and Office of the Vice President for Information Technology

This analysis was performed on the [zcluster](https://wiki.gacrc.uga.edu/wiki/Systems) at the University of Georgia

### System information

```
$ uname -a

Linux zcluster.rcc.uga.edu 2.6.18-398.el5 #1 SMP Tue Aug 12 06:26:17 EDT 2014 x86_64 x86_64 x86_64 GNU/Linux

$ lsb_release --all

LSB Version:	:core-4.0-amd64:core-4.0-ia32:core-4.0-noarch:graphics-4.0-amd64:graphics-4.0-ia32:graphics-4.0-noarch:printing-4.0-amd64:printing-4.0-ia32:printing-4.0-noarch
Distributor ID:	RedHatEnterpriseServer
Description:	Red Hat Enterprise Linux Server release 5.11 (Tikanga)
Release:	5.11
Codename:	Tikanga
```

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
	- this is required only to check the read quality before analysis



### Prepare environment
1. add fastq files to `ESSENTIAL/FASTQ`
	- follow naming conventions specified in `ESSENTIAL/FASTQ/naming_conventions.txt`
2. specify reference parent and divergent parent filenames in the `ESSENTIAL/FASTQ/parents.txt`
3. edit `ESSENTIAL/SCRIPTS/020.compare_data_to_genomes.sh` to run proper mapping script
	- this git repo contains scripts for mapping paired end and single end reads
	- if running paired end script, remove the `#` pound sign comment before `sh $FILES/021.paired_end_map_to_OR.sh`
		- also, edit the min and max fragment size parameters in the `ESSENTIAL/SCRIPTS/021.paired_end_map_to_OR.sh` file for the bowtie2 command (line 65)
	- if running single end script, remove the `#` pound sign comment before `sh $FILES/021.single_end_map_to_OR.sh`=

### steps to run
1. run `010.check_data.sh`
	- check output data in `FASTQC_OUT`, and clean data as needed
		- i.e. trim reads, k-mer balance with [Lighter](https://github.com/mourisl/Lighter), etc.
2. run `020.compare_data_to_genomes.sh` & `030.determine_snp_profiles_for_samples.sh`
3. graph output data from snp mapping located in `SNP_MAPPING/PARSED_SNP_INFO`
	- make a scatter plot of RATIO ~ KB for each CONTIG, locate regions of interest
4. list your regions of interest in .bed format files in the `ESSENTIAL/FILTER_SITES` directory
	- follow guidelines listed in the readme.txt file in that folder
5. run `040.check_for_synonymy.sh`
	- obtain output in `GFF_OVERLAP` and `GFF_OVERLAP/TRANSLATE_CDS` folders
	- `GFF_OVERLAP/*.all` lists all GFF features that high quality snps overlap
		- format: {full GFF entry}{full vcf entry for snp}
	- `GFF_OVERLAP/*.snps_not_in_genes` lists all snps that fall in euchromatin regions but do not overlap any GFF features
		- these may hit promotors or other factors outside the gene body, and may be of interest if a mutation is not found in the GFF features
	- `GFF_OVERLAP/TRANSLATE_CDS/*.translated_CDS` outputs the translated AA sequence for all snps falling in coding sequences, and shows if the snp
		produces a non-synonymous output
