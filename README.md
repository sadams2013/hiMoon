# hiMoon: Determine Named Haplotypes

![himoon_hap](https://github.com/Ariel-Precision-Medicine/himoon_hap/workflows/himoon_hap/badge.svg)

## Disclaimer and Status

Haplotype calling and structural variant estimation is highly dependent on the quality of data provided. 
The method is also under development, and will change. 
This can cause results to change as it undergoes further testing and development. 
In short: this will produce results, but those results might not be right. 

Intended for evaluation and research use only.
User assumes all risk associated with data generated from this platform. 

## Description

hiMoon is a diversly applicable software tool that estimates named haplotypes (e.g. star alleles). 
The limits of the tool are driven by the needs of the end user. 
It uses argued translation tables in the format used by PharmVar for the determination of haplotypes. 

## Installation

I recommend standalone installation with Anaconda. 

```
conda env create -q -n himoon -f test_env.yml
conda activate himoon
python setup.py install
```

## Operation of hiMoon

### Configuration

Configuration information for hiMoon is specified in a file called config.ini, which must be located in the directory from which hiMoon is being called. 
If one is not present, hiMoon will create one with default values. 
The config.ini file in this repository is commented with all possible fields to document structure. 

### Running from the command line

After installation, hiMoon can be invoked by simply calling: 

```
hiMoon -h

usage: hiMoon [-h] [-f VCFFILE] [-t TRANSLATION_TABLES] [-o OUTPUT_DIRECTORY]
              [-c CONFIG_FILE] [-i] [-s SAMPLE]

Match haplotypes, return raw data and/or reports.

optional arguments:
  -h, --help            show this help message and exit
  -f VCFFILE, --vcffile VCFFILE
                        path/to/vcf file
  -t TRANSLATION_TABLES, --translation-tables TRANSLATION_TABLES
                        Directory with translation tables or a single
                        translation table file
  -o OUTPUT_DIRECTORY, --output-directory OUTPUT_DIRECTORY
                        Directory for Output Files. If not specified, will
                        output calls to stdout.
  -c CONFIG_FILE, --config-file CONFIG_FILE
                        path to config file.
  -i, --loglevel-info   Use more verbose logging output (useful for
                        debugging).
  -s SAMPLE, --sample SAMPLE
                        Single sample from multisample ID (if not specified,
                        will do all)
```

You must provide a VCF file that is indexed. 


### Method for Haplotype Determination

hiMoon operates under the assumption of orthogonal haplotype assignment. 
This approach is relatively conservative and is most appropriate for instances where discovery of novel alleles is not desired. 
Haplotypes are estimated in a stepwise process. 
1. For each haplotype in a translation table, iterate over all associated markers and flag if they match the haplotype. 
2. Optimize such that two haplotypes are matched with the maximum amount of found variants on the translation table. 

## Contributing

I welcome any contributions, corrections, comments, criticisms. 
