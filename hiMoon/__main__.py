import argparse
import os
import glob
import sys
import csv

from .subject import Subject
from .gene import AbstractGene
from .vcf import VarFile, write_variant_file

from . import LOGGING, CONFIG, set_config, set_logging_info

def get_vcf_genes(args) -> ([AbstractGene], VarFile):
    """
    Prep VCF and gene objects

    Args:
        args ([type]): args

    Returns:
        Tuple
    """
    vcf = VarFile(args["vcffile"], args["sample"])
    genes = []
    if args["translation_tables"][-3:] == "tsv":
        genes.append(AbstractGene(os.path.abspath(args["translation_tables"]), vcf))
    else:
        for translation_table in glob.glob(args["translation_tables"] + "/*.tsv"):
            genes.append(AbstractGene(os.path.abspath(translation_table), vcf))
    return vcf, genes

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Match haplotypes, return raw data and/or reports.", prog="hiMoon")
    parser.add_argument("vcffile", help="path/to/vcf file")
    parser.add_argument("-t", "--translation-tables",
                        help="Directory with translation tables or a single translation table file", 
                        default=None)
    parser.add_argument("-o", "--output-directory",
                        default="./",
                        help="Directory for Output Files.")
    parser.add_argument("-c", "--config-file",
                        default=None,
                        help="path to config file.")
    parser.add_argument("-i", "--loglevel-info",
                        help="Use more verbose logging output (useful for debugging).",
                        action="store_true")
    parser.add_argument("-s", "--sample",
                        help="Single sample from multisample ID (if not specified, will do all)",
                        default=None)
    
    args = vars(parser.parse_args())
    if args["translation_tables"] is None:
        print("You must provide a translation table or a directory with translation tables.")
        sys.exit(1)
    if args["loglevel_info"]:
        set_logging_info()
    if args["config_file"]:
        set_config(args["config_file"])
    vcf, genes = get_vcf_genes(args)
    subjects = [Subject(prefix = sub_id, genes = genes) for sub_id in vcf.samples]
    out_dir = args["output_directory"]
    write_variant_file(out_dir, subjects, args["vcffile"].split("/")[-1].strip(".vcf.gz").strip(".bcf"), genes)

if __name__ == "__main__": 
    main()
    