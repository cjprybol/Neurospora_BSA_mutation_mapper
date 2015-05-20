# Neurospora BSA mutation mapper
---

## About
This project was created as an analysis pipeline to locate causitive suppressor mutations in the model organism Neurospora crassa.

### setup project
To use this pipeline on your own data, or to recreate the findings out our study,
start by cloning this repository into a directory called `ESSENTIAL`

```
mkdir {desired working directory}
cd {desired working directory}
git clone https://github.com/cprybol/DIM5_suppressor_mapping.git ESSENTIAL

```

### required software (all must be in user path)
- Perl 5
- Python2
	- must be the default system python (required by Bowtie for mapping reads to genome)
- Python3
	- recommended installation method is to simply get [Anaconda](http://continuum.io/downloads#34)
		- restart bash session after initial installation
	- **packages**	
		- [pandas](http://pandas.pydata.org/)
			- included in Anaconda
		- [BioPython](http://biopython.org/wiki/Main_Page)
			- not included in Anaconda
			- acquire by running `conda install biopython`
- R
	- version 3+
	- **packages**
		- [ggplot2](http://ggplot2.org/)
		- [gridExtra](http://cran.r-project.org/web/packages/gridExtra/index.html)
- [Samtools](http://www.htslib.org/)
- [bcftools](http://www.htslib.org/)
- [htslib](http://www.htslib.org/)
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [Bedtools](http://bedtools.readthedocs.org/en/latest/)
- [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [Lighter](https://github.com/mourisl/Lighter)
- [BreakDancer](http://gmt.genome.wustl.edu/packages/breakdancer/)
- [Pindel](http://gmt.genome.wustl.edu/packages/pindel/)


### Prepare environment
1. add fastq files to the directory `{path to working directory}/ESSENTIAL/FASTQ`
	- follow naming conventions specified in `{path to working directory}/ESSENTIAL/FASTQ/naming_conventions.txt`
	- fastq files must be gzip-ed
2. specify reference parent and divergent parent filenames, read-type (single- or paired-end), # of available cores, and min and max library fragment sizes in the `{path to working directory}/ESSENTIAL/config.txt` file

### steps to run
1. run `{path to working directory}/ESSENTIAL/SCRIPTS/001.master.sh`
2. evaluate scatterplots in `{path to working directory}/SNP_MAPPING/PARSED_SNP_INFO/GRAPHS` and locate regions that satisfy desired similarity thresholds
3. list your regions of interest in .bed format files in the `{path to working directory}/ESSENTIAL/FILTER_SITES` directory
	- follow guidelines listed in the `{path to working directory}/ESSENTIAL/FILTER_SITES/readme.txt` file
4. run `{path to working directory}/ESSENTIAL/SCRIPTS/002.master.sh`
	- obtain output in `{path to working directory}/GFF_OVERLAP` and `{path to working directory}/GFF_OVERLAP/TRANSLATE_CDS` folders
	- `{path to working directory}/GFF_OVERLAP/*.all` lists all GFF features that high quality snps overlap
		- format: {full GFF entry}{full vcf entry for snp}
	- `{path to working directory}/GFF_OVERLAP/*.snps_not_in_genes` lists all snps that fall in euchromatin regions but do not overlap any GFF features
		- these may hit promotors or other factors outside the gene body, and may be of interest if a mutation is not found in the GFF features
	- `{path to working directory}/GFF_OVERLAP/TRANSLATE_CDS/*.translated_CDS` outputs the translated AA sequence for all snps falling in coding sequences, and shows if the snp produces a non-synonymous output
	- `{path to working directory}/SV_DETECTION` contains BreakDancer and Pindel output to detect possible structural variants
