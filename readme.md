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


### steps to run

All necessary elements for running the analysis on this data are located in the ESSENTIAL directory
In order to reproduce the full data set and analysis, run the scripts found in ESSENTIAL/SCRIPTS in numeric order
NOTE: the 'BASE' directory will need to be changed to reflect the directory location in which the ESSENTIAL folder is placed
