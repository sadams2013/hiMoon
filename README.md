# hiMoon: Determine Pharmacogenomics Alleles and Structural Variation

[![Build Status](https://travis-ci.org/sadams2013/hiMoon.svg?branch=master)](https://travis-ci.org/sadams2013/hiMoon)

## Disclaimer and Status

**This is a project in ongoing active development that is currently in pre-alpha.**

Haplotype calling and structural variant estimation is highly dependent on the quality of data provided. 
The method is also under development, and will change. 
This can cause results to change as it undergoes further testing and development. 
In short: this will produce results, but those results might not be right. 

Intended for evaluation and research use only.
User assumes all risk associated with data generated from this platform. 

## Why make this?

- Pharmacogenomics platforms are diverse (e.g. array, WES, WGS) and each platform is variable. 
- Applications have different needs that require tuning of software tools. 
- Efforts by existing authorities (CPIC, PharmVar, PharmGKB) to standardize definitions are sporadically implemented. 
- An existing open source solution for pharmacogenomics haplotyping and structural variation determination does not exist. 
    - Note: PharmCAT relies on Astrolabe. which is closed source and not available without cost to non-academic entities.

## Description

hiMoon is a diversly applicable software tool that estimates pharmacogenomics haplotypes (e.g. star alleles) and copy number variation. 
The limits of the tool are driven by the needs of the end user. 
It uses argues translation tables in the format used by PharmVar for the determination of haplotypes. 
Regions of interest for structural variation are provided in a configuration file. 

## Installation

hiMoon has been tested on MacOS and Linux. 
It will probably work on Windows, but might need some tweaks. 

### Dependencies

- numba
- pysam
- numpy

### Install from source

It is highly recommended that you install hiMoon in an isolated Python environment (e.g. virtualenv). 

git clone https://github.com/sadams2013/hiMoon.git
pip install -r requirements.txt
python setup.py install

## Operation of hiMoon

### Configuration

Configuration information for hiMoon is specified in a file called config.ini, which must be located in the directory from which hiMoon is being called. 
If one is not present, hiMoon will create one with default values. 
The config.ini file in this repository is commented with all possible fields to document structure. 

### Running from the command line

After installation, hiMoon can be invoked by simply calling: 

```
hiMoon -h

usage: hiMoon [-h] [-f VCFFILE] [-b BAMFILE] [-t TRANSLATION_TABLES]
              [-c COPY_NUMBER_CONTROL] [-i] [-o OUTPUT_DIRECTORY]
              [-p OUTPUT_PREFIX]

Match haplotypes, return raw data and/or reports.

optional arguments:
  -h, --help            show this help message and exit
  -f VCFFILE, --vcffile VCFFILE
                        path/to/vcf file
  -b BAMFILE, --bamfile BAMFILE
                        path/to/bam file (also needs to be indexed)
  -t TRANSLATION_TABLES, --translation-tables TRANSLATION_TABLES
                        Directory with translation tables or a single
                        translation table file
  -c COPY_NUMBER_CONTROL, --copy-number-control COPY_NUMBER_CONTROL
                        BAM file with reference copy number.
  -i, --loglevel-info   Use more verbose logging output (useful for
                        debugging).
  -o OUTPUT_DIRECTORY, --output-directory OUTPUT_DIRECTORY
                        Directory for Output Files. If not specified, will
                        output calls to stdout.
  -p OUTPUT_PREFIX, --output-prefix OUTPUT_PREFIX
                        prefix for output files.
```

You must provide either a sample VCF file or a sample BAM file. 
You can also provide both a VCF and a BAM file. 
If both files are provided, hiMoon will use its variant caller if a for genes specified in config.ini. 
If only a VCF is provided, hiMoon will use variants from that for all translation tables provided. 
If only a BAM is provided, hiMoon will use its allele caller for all translation tables. 
*Note: most variant callers (e.g. GATK) will do a better job than hiMoon on a large scale.* 
*hiMoon's variant caller is only really appropriate if you want finer control over specific regions.*

If you want to estimate copy number variation, you must provide a BAM file for the sample. 
If you only provide a BAM file for the sample, then hiMoon will attempt to use 'Quick CNV.'
This is only intended to serve as an estimator, and is a relatively crude approach. 
'Quick CNV' relies on a calibrator region in the sample file. 
We recommend a region in *CYP2D8P* with consistent coverage, which is used to calculate a normalized coverage ratio. 

To use the maximum penalized likelihood estimate (MPLE), you must also provide a control BAM file. 
This file should be sequenced with the same methodology/chemistry as the unknown sample. 
The control file must also be known to not have any structural variation. 

### Method for Haplotype Determination

hiMoon operates under the assumption of orthogonal haplotype assignment. 
This approach is relatively conservative and is most appropriate for instances where discovery of novel alleles is not desired. 
hiMoon also returns a boolean estimate of whether or not a sample appears to have a novel haplotype. 
It should be noted that this refers to the incidental finding, and does not reflect an exhaustive search for novel haplotypes. 
It is instead intended to suggest that users investigate base-pair level data rather than rely on the called haplotypes. 
Haplotypes are estimated in a stepwise process. 
1. For each haplotype in a translation table, iterate over all associated markers and flag if they match the haplotype. 
2. Iterate over all haplotypes with 100% identity with the sample. 

This method is robust in instances of complex variation (e.g. CYP2D6). 
It is also relatively robust to gaps in data and will note ambiguities in coverage. 


### Method for Structural Variation Determination

Structural variation is determined for NGS data (bam file required for sample and a control) with the maximum penalized likelihood estimate method using depth of coverage adapted from the method described by xxx and colleagues. 
However, rather than looking at all regions for capture based sequencing, the user provides target regions that span large regions of coverage. 
hiMoon iteratively and recursively determines the ideal starting point (maximized likelihood) for each candidate break point defined by the start of an individual read. 
Likelihood maximized regions with a coverage ratio > 1.4 are reported as copy number gain, and < 0.6 are reported as copy number loss. 

## Contributing

I welcome any contributions, corrections, comments, criticisms. 


## TODO

- More expansive documentation.
- Large scale test of copy number estimate in WGS, WES, and targeted sequencing data. 
- Ability to call CNV loss and gain directly in results. 