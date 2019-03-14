import argparse
import os
import glob
import sys
import csv

from .subject import Subject
from .gene import Gene

from . import logging


def main():
    parser = argparse.ArgumentParser(
        description="Match haplotypes, return raw data and/or reports.", prog="hiMoon")
    parser.add_argument("-f", "--vcffile",
                        help="path/to/vcf file", default="")
    parser.add_argument("-b", "--bamfile",
                        help="path/to/bam file (also needs to be indexed)", default="")
    parser.add_argument("-t", "--translation-tables",
                        help="Directory with translation tables or a single translation table file",
                        default=None)
    parser.add_argument("-c", "--copy-number-control",
                        default=None,
                        help="BAM file with reference copy number.")
    parser.add_argument("-i", "--loglevel-info",
                        help="Use more verbose logging output (useful for debugging).",
                        action="store_true")
    parser.add_argument("-o", "--output-directory",
                        default="./",
                        help="Directory for Output Files. \
                        If not specified, will output calls to stdout.")
    parser.add_argument("-p", "--output-prefix",
                        default="out",
                        help="prefix for output files.")
                        
    args = vars(parser.parse_args())

    # The log can get pretty verbose if translation tables are very long. 
    if args["loglevel_info"]:
        logging.getLogger().setLevel(logging.INFO)

    if args["vcffile"][-6:] == "vcf.gz" or args["bamfile"]:
        haplotyper(args)
    else:
        print("You must specify at least one input file (vcf or bam)")
        sys.exit(1)


def haplotyper(args: dict):
    """
    """
    genes = []
    if args["bamfile"]:
        bam_file = os.path.abspath(args["bamfile"])
    else:
        bam_file = None
    if args["vcffile"]:
        vcf_file = os.path.abspath(args["vcffile"])
    else:
        vcf_file = None
    if args["copy_number_control"]:
        copy_number_control = os.path.abspath(args["copy_number_control"])
    else:
        copy_number_control = None

    # Test if single translation table argued, or directory of files. 
    if args["translation_tables"][-3:] == "tsv":
        genes.append(Gene(os.path.abspath(args["translation_tables"])))
    else:
        for translation_table in glob.glob(args["translation_tables"] + "/*.tsv"):
            genes.append(Gene(os.path.abspath(translation_table)))
    subj = Subject(prefix = args["output_prefix"], genes = genes, vcf_file=vcf_file, alignment_file=bam_file, copy_number_control=copy_number_control)
    if args["output_directory"]:
        out_dir = args["output_directory"]
        write_report(out_dir, subj, args["output_prefix"])
        write_cnv_report(out_dir, subj, args["output_prefix"])
    else:
        output_to_stdout(subj)



def output_to_stdout(subject: Subject):
    """
    Dump output to stdout - allows piping or looping over several files
    """
    for gene, data in subject.called_haplotypes.items():
        print(f"{gene}\t{data[1]}\t"
              f"{'/'.join([str(hap) for hap in data[0]])}\t"
              f"{'|'.join([','.join([str(hap) for hap in dip]) for dip in data[2]])}\t"
              f"{data[3]}")

def write_report(directory: str, subject: Subject, prefix):
    """
    """
    rows = []
    for gene, data in subject.called_haplotypes.items():
        rows.append(
            {
                "PREFIX": subject.prefix,
                "GENE": gene,
                "VERSION": data[0],
                "CALLED": "/".join([str(hap) for hap in data[1]]),
                "POSSIBLE": "|".join([",".join([str(hap) for hap in dip]) for dip in data[2]]),
                "NOVEL_HAPLOTYPE_DETECTED": data[3]
            }
        )
    with open(directory+f"/{prefix}.haplotypes.tsv", "w") as outfile:
        fieldnames = ["PREFIX", "GENE", "VERSION", "CALLED", "POSSIBLE", "NOVEL_HAPLOTYPE_DETECTED"]
        outwriter = csv.DictWriter(outfile, delimiter="\t", fieldnames=fieldnames)
        outwriter.writeheader()
        outwriter.writerows(rows)

def write_cnv_report(directory: str, subject: Subject, prefix):
    with open(directory + f"/{prefix}.cnv.tsv", "w") as cnvout:
        fieldnames = ["PREFIX", "REGION", "CHROMOSOME", "SPAN_START", "SPAN_STOP", "LENGTH", "SPAN_RATIO"]
        outwriter = csv.DictWriter(cnvout, delimiter="\t", fieldnames=fieldnames)
        outwriter.writeheader()
        outwriter.writerows(subject.cnv_regions)


if __name__ == "__main__": 
    main()