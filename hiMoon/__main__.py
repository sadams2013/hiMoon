import argparse
import os
import glob
import sys
import csv

from .subject import Subject
from .gene import Gene
from .config import ConfigData
from .vcf import VarFile

from . import logging


def main():
    parser = argparse.ArgumentParser(
        description="Match haplotypes, return raw data and/or reports.", prog="hiMoon")
    parser.add_argument("-f", "--vcffile",
                        help="path/to/vcf file", default="")
    parser.add_argument("-t", "--translation-tables",
                        help="Directory with translation tables or a single translation table file",
                        default=None)
    parser.add_argument("-o", "--output-directory",
                        default="./",
                        help="Directory for Output Files. \
                        If not specified, will output calls to stdout.")
    parser.add_argument("-c", "--config-file",
                        default="config.ini",
                        help="path to config file.")
    parser.add_argument("-i", "--loglevel-info",
                        help="Use more verbose logging output (useful for debugging).",
                        action="store_true")
    parser.add_argument("-s", "--sample",
                        help="Single sample from multisample ID (if not specified, will do all)",
                        default=None)
    
    args = vars(parser.parse_args())

    # The log can get pretty verbose if translation tables are very long. 
    if args["loglevel_info"]:
        logging.getLogger().setLevel(logging.INFO)

    haplotyper(args)



def haplotyper(args: dict):
    """
    """
    config = ConfigData(args["config_file"])
    genes = []
    subjects = []
    vcf = VarFile(args["vcffile"], args["sample"])

    # Test if single translation table argued, or directory of files. 
    if args["translation_tables"][-3:] == "tsv":
        genes.append(Gene(os.path.abspath(args["translation_tables"]), config, vcf))
    else:
        for translation_table in glob.glob(args["translation_tables"] + "/*.tsv"):
            genes.append(Gene(os.path.abspath(translation_table), config))
    contigs = [g.chromosome for g in genes]
    subjects = [Subject(prefix = sub_id, genes = genes) for sub_id in vcf.samples]
    out_dir = args["output_directory"]
    write_report(out_dir, subjects, args["vcffile"].split("/")[-1])


def write_report(directory: str, subjects: [Subject], prefix):
    """
    """
    rows = []
    for subj in subjects:
        for gene, data in subj.called_haplotypes.items():
            rows.append(
                {
                    "PREFIX": subj.prefix,
                    "GENE": gene,
                    "VERSION": data[0],
                    "CALLED": "/".join([str(hap) for hap in data[3]]),
                    "MARKERS": ", ".join(data[2]),
                }
            )
    with open(directory+f"/{prefix}.haplotypes.tsv", "w") as outfile:
        fieldnames = ["PREFIX", "GENE", "VERSION", "CALLED", "MARKERS"]
        outwriter = csv.DictWriter(outfile, delimiter="\t", fieldnames=fieldnames)
        outwriter.writeheader()
        outwriter.writerows(rows)


if __name__ == "__main__": 
    main()