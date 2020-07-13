# hiMoon: Determine Named Haplotypes from VCF

![hiMoon](https://github.com/sadams2013/hiMoon/workflows/himoon_hap/badge.svg)
![PyPI](https://github.com/sadams2013/hiMoon/workflows/Upload%20to%20PyPI/badge.svg)

## Read before using

hiMoon is provided such that you might find it useful. 
Users should be aware that named haplotype calling is dependent on the breadth and quality of data provided. 
If you use this tool, you need to understand how it works and inherent limitations of the methods employed. 

### Main Rules for Use

1. Know your input data! hiMoon will work with what you give it, so Garbage In Garbage Out applies. 
2. Function will be erratic until a 1.0 release (i.e. beta only), this is especially important for CYP2D6. 
3. Validate with known samples + edge cases before using in any meaningful application. 

**hiMoon is intended for research use only.**

## Description

hiMoon is a diversly applicable software tool that estimates named haplotypes (e.g. star alleles). 
It uses argued translation tables in the format used by PharmVar for the determination of haplotypes. 

## Useage 

### Installation

Install with pip: 

```pip install hiMoon```

Or clone this repository and run: 

```python setup.py install```

### Configuration

hiMoon configuration is managed through config.py, which provides several default configuration parameters that coveres most simple use-cases. 
However, custom configuration can be provided by defining a file called config.ini. 
This file can define contig accessions, VCF parsing parameters, and IUPAC codes. 
While it is not expected that these will need to be manually defined/overridden in most cases, there are certainly instances where users might want to tweak them. 

To create a base configuration file based on the default parameters, you can call ```hiMoon -c default```, which will write a file called 'himoon_config.ini' that you can modify as needed. 

### hiMoon CLI

Working with the command line interface is fairly straightforward. 
You can start with: 

```
usage: hiMoon [-h] [-t TRANSLATION_TABLES] [-o OUTPUT_DIRECTORY]
              [-c CONFIG_FILE] [-i] [-s SAMPLE] [-S SOLVER]
              [-M ALLOWED_NO_MATCH]
              [vcf_file]

Match haplotypes, return raw data and/or reports.

positional arguments:
  vcf_file              path/to/vcf file

optional arguments:
  -h, --help            show this help message and exit
  -t TRANSLATION_TABLES, --translation-tables TRANSLATION_TABLES
                        Directory with translation tables or a single
                        translation table file
  -o OUTPUT_DIRECTORY, --output-directory OUTPUT_DIRECTORY
                        Directory for Output Files.
  -c CONFIG_FILE, --config-file CONFIG_FILE
                        path to config file. -c default will write a config
                        file that can be modified for this.
  -i, --loglevel-info   Use more verbose logging output (useful for
                        debugging).
  -s SAMPLE, --sample SAMPLE
                        Single sample from multisample ID (if not specified,
                        will do all)
  -S SOLVER, --solver SOLVER
                        Solver to use (GLPK or CBC), default = CBC
```

You must provide a compressed (.vcf.gz, .bcf) and indexed (.tbi, .csi) VCF file. 
Ideally, your VCF should be a genomic vcf (gVCF, g.vcf), meaning that all callable positions are defined.
hiMoon assumes that any missing positions in the VCF file are "no-call."
Therefore, discretely defined genotypes will increase the specificity of named-haplotype calls. 

You must also provide a translation table (-t, --translation-tables). 
This argument can be to a specific file, or to a directory containing at least one translation table. 
The format for these tables is the same format that can be exported from PharmVAR. 
While several genes are available from PharmVAR, users will also need to manually produce translation tables for their genes of interest. 
See the "Translation Table Format" section for details on the format. 

If you provide a multi-sample VCF, hiMoon will automatically perform haplotype matching for each sample in your VCF. 
To select a single sample from the VCF, use ```-s SAMPLE_ID``` to select your specific sample. 

#### Output

By default, hiMoon produces a valid VCF v4.3 that contains per sample haplotype calls. 
The REF column is populated with the reference haplotype (e.g. *1), and the alts correspond to named haplotypes. 
For multi-sample VCF files and complex genes (e.g. CYP2D6), this can produce a substantial ALT column. 
A "HC" score is provided for each sample that corresponds to 1 / {total number of possible genotypes}. 
**Anytime the HC score is < 1 should trigger further investigation!**

hiMoon will alo create a TSV file that has one sample + gene call per line. 
If multiple possible haplotype combinations are found, each call will be on a separate line. 

### API

hiMoon is also exposed through a simple API. 

```python
from hiMoon import himoon

get_haps_from_vcf(translation_table_path: str, vcf_file_path: str, sample_id: str, config_path = None) -> tuple

get_haps_from_variants(translation_table_path: str, vcf_data: str, sample_id: str, config_path = None) -> tuple
```

## Issues

Please use issues in this repo to report issues or bugs. 

## Contributing

Contributions through pull requests are welcome. 
Also criticisms, suggestions, etc... through issues are also fine with me. 

# Appendices

## Translation Table Format

All translation table files are formatted per the exportable haplotype definitions from PharmVAR (.tsv) files. 
A .tsv haplotype definition file is required for hiMoon. 
Optionally, users can define a separate .cnv file that must be otherwise named exactly as the main .tsv file. 

Prospective users are encouraged to look at example tables in hiMoon/tests/test_files/translation_tables for full examples of table formatting.
The parsing methods for these tables are in the AbstractGene class of gene.py, and can provide further insight into how one can format these tables. 

### Haplotype Definition Files (.tsv)

This file is read into hiMoon as a whitespace-delimited file. 
An example of formatting is shown below: 

```
#version=pharmvar-4.1.6
Haplotype Name	Gene	rsID	ReferenceSequence	Variant Start	Variant Stop	Reference Allele	Variant Allele	Type
CYP2D6*1	CYP2D6		REFERENCE	.				
CYP2D6*1.001	CYP2D6		REFERENCE	.				
CYP2D6*1.002	CYP2D6	rs28371732	NC_000022.11	42126963	42126963	C	T	substitution
CYP2D6*1.003	CYP2D6	rs150163869	NC_000022.11	42128813	42128813	G	A	substitution
...
CYP2D6*139	CYP2D6	rs1058172	NC_000022.11	42127526	42127526	C	T	substitution
CYP2D6*139.001	CYP2D6	rs28371713	NC_000022.11	42128793	42128793	A	G	substitution
CYP2D6*139.001	CYP2D6	rs374672076	NC_000022.11	42131189	42131189	C	G	substitution
CYP2D6*139.001	CYP2D6	rs113889384	NC_000022.11	42128741	42128741	G	A	substitution
CYP2D6*139.001	CYP2D6	rs112568578	NC_000022.11	42128711	42128711	C	G	substitution
CYP2D6*139.001	CYP2D6	rs111564371	NC_000022.11	42128706	42128706	T	C	substitution
CYP2D6*139.001	CYP2D6	rs4987144	NC_000022.11	42127001	42127001	G	A	substitution
CYP2D6*139.001	CYP2D6	rs28371726	NC_000022.11	42127537	42127537	A	G	substitution
CYP2D6*139.001	CYP2D6	rs1058172	NC_000022.11	42127526	42127526	C	T	substitution
```

Reading begins at the third line. 
The version information is also extracted and stored separately. 
Multiple REFERENCE rows can be defined, but only the first appearance is used. 
In this case, if one or zero named haplotypes are matched, hiMoon will report *1 since it is defined as the reference. 
If no reference is defined, hiMoon will default to "REF" for any non matching alleles.

Main + sub-alleles are denoted with a "." following the main allele. 
Main alleles with and without sub-alleles defined can be provided in the same translation table file, or users may elect to modify the table to only contain the main allele or main + sub-alleles. 
Considering that main + sub-alleles are more specific, it is likely prudent to use those with data that is very broad (e.g. g.vcf for a WGS), whereas only certain main alleles might be provided if data are sparse (e.g. targeted sequencing, non g.VCF, targeted array).

### Structural Variant Definition File (.cnv)

While users can use the .tsv file to define structural variants that are to be called as named haplotypes, it is likely preferable to define those in a separate .cnv file. 
Formatting for this table is the exact same as the .tsv haplotype definition table, but there are distinctions in how it is parsed. 
These variants use the "CNV" type. 
This is highly consequential in that the optimizer is forced to use any found CNV variants. 
Contact the maintainer of hiMoon if you want to incorporate SVs that are weighted the same as other variants. 

An example is provided below: 

```
#version=cyp2d6_GRCh37_1000_genomes
Haplotype Name	Gene	rsID	ReferenceSequence	Variant Start	Variant Stop	Reference Allele	Variant Allele	Type
CYP2D6*5.001	CYP2D6	CYP2D6	NC_000022.10	42523949	42533891	T	<CN0>	CNV
CYP2D6*2_x2	CYP2D6	CYP2D6	NC_000022.10	42523949	42533891	T	<CN2>	CNV
CYP2D6*4_x2	CYP2D6	CYP2D6	NC_000022.10	42523949	42533891	T	<CN2>	CNV
```

In this example, we define two copy-number gain alleles and a copy-loss variant. 
CYP2D6*5 defines a whole gene deletion and is quite straightforward to define. 
However, *2xN and *4xN are defined as copy gain + variations. 
Therefore, calling them with this method requires that we define all main alleles + suballeles + CNV. 
The CNV file is parsed such that (in this case), *5.001 is simply appended to the haplotype table. 
*2xN and *4xN are expanded to all *2 and *4 sub-alleles present in the haplotype table. 
Since *4 and *2 are represented by hundreds of permutations of variants across tens of sub-alleles, the full definition of these alleles + copy gain is quite cumbersome. 
This method simplifies that substantially. 

Note that the decision to include *5 in the haplotype file vs. the CNV file is arbitrary, but I recommend it to keep SVs organized and to avoid unnecessary modifications to the tsv file. 

#### How to define SVs

hiMoon is designed to work with the data provided, but for better or for worse, it will do what you ask it to do. 
Methods for SV determination and reporting vary, and this can have major implications on how you define SVs in the cnv file. 

For example: the hiMoon test files are based on data reported in Phase 3 of the 1000 genomes project. 
One set of tests looks for correct named haplotypes (based on Getrm consensus) in GRCh38 low-depth VCF calls with no SVs defined. 
The other looks at a combined SNPs, InDels, and SVs from the GRCh37 calls with a defined .cnv file (see example above). 
The breakpoints in that file are defined based on manual inspection of the structural variant VCF, wherein we found that the CNV for *CYP2D6* is defined as starting as position 42523949. 
We also noted that this is defined for CN0, CN2, CN3, and CN4. 
Therefore, our .cnv table MUST match these parameters. 

Different methods for investigating CNVs might have placed the breakpoint before or after this defined point, so we would need to verify that and modify if this table were used for a different dataset. 
Furthermore, copy gain/loss could be defined by a different alt allele (e.g. "DEL", "DUP"). 
This is especially important for more complex alleles, like CYP2D6/7 hybrid alleles (e.g. \*4+\*68), which may or may not be determinable depending on the upstream method for finding SVs.
Of note, the common NA12878 sample in 1000 genomes is \*4+\*68/\*3, but would be determined using this dataset to be \*4x2/\*3. 
This is a limitation of SV calling methods and not a limitation of hiMoon. 
