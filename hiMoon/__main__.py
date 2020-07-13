#    Copyright 2020 Solomon M. Adams

#    Licensed under the Apache License, Version 2.0 (the "License");
#    you may not use this file except in compliance with the License.
#    You may obtain a copy of the License at

#        http://www.apache.org/licenses/LICENSE-2.0

#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#    See the License for the specific language governing permissions and
#    limitations under the License.

import argparse
import os
import glob
import sys
import csv

from .subject import Subject
from .gene import AbstractGene
from .vcf import VarFile, write_variant_file, write_flat_file

from . import LOGGING, CONFIG, set_config, set_logging_info

def get_vcf_genes(args) -> ([AbstractGene], VarFile):
    """
    Prep VCF and gene objects

    Args:
        args ([type]): args

    Returns:
        Tuple
    """
    vcf = VarFile(args["vcf_file"], args["sample"])
    genes = []
    solver = args["solver"]
    if args["translation_tables"][-3:] == "tsv":
        genes.append(AbstractGene(os.path.abspath(args["translation_tables"]), vcf, solver = solver))
    else:
        for translation_table in glob.glob(args["translation_tables"] + "/*.tsv"):
            genes.append(AbstractGene(os.path.abspath(translation_table), vcf, solver = solver))
    return vcf, genes

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Match haplotypes, return raw data and/or reports.", prog="hiMoon")
    parser.add_argument("vcf_file", help="path/to/vcf file", nargs="?")
    parser.add_argument("-t", "--translation-tables",
                        help="Directory with translation tables or a single translation table file", 
                        default=None)
    parser.add_argument("-o", "--output-directory",
                        default="./",
                        help="Directory for Output Files.")
    parser.add_argument("-c", "--config-file",
                        default=None,
                        help="path to config file. -c default will write a config file that can be modified for this.")
    parser.add_argument("-i", "--loglevel-info",
                        help="Use more verbose logging output (useful for debugging).",
                        action="store_true")
    parser.add_argument("-s", "--sample",
                        help="Single sample from multisample ID (if not specified, will do all)",
                        default=None)
    parser.add_argument("-S", "--solver",
                        help="Solver to use (GLPK or CBC), default = CBC",
                        default="CBC")
    
    args = vars(parser.parse_args())
    if args["config_file"] ==  "default":
        LOGGING.warning("Writing default config file and exiting")
        CONFIG.write_config("himoon_config.ini")
        sys.exit(0)

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
    write_variant_file(out_dir, subjects, args["vcf_file"].split("/")[-1].strip(".vcf.gz").strip(".bcf"), genes)
    write_flat_file(out_dir, subjects, args["vcf_file"].split("/")[-1].strip(".vcf.gz").strip(".bcf"))

if __name__ == "__main__": 
    main()
    