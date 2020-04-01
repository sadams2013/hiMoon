"""
Defines the class Gene and Haplotype
Gene contains methods useful for the gene
Haplotype contains methods for gene specific haplotypes
"""

import pandas as pd
import sys

from . import logging
from .vcf import VarFile

class Gene:
    """
    Describes a gene object and provides clean library data.
    This contains pertinent information about the gene
    and its haplotypes
    """

    def __init__(self, translation_table, config, vcf: VarFile):
        self.matched = False
        self.haplotypes = {}
        self.gene = None
        self.accession = None
        self.ref = "Ref"
        self.config = config
        with open(translation_table, 'rt') as trans_file:
            self.version = trans_file.readline().strip("#version=\n\t")
        self.translation_table = pd.read_csv(
            translation_table, 
            skipinitialspace = True,
            delimiter = "\t", 
            skiprows = 1,
            na_values= {4: ".", 5: "."},
            dtype = {4: pd.Int64Dtype(), 5: pd.Int64Dtype()})
        self.gene = self.translation_table.iloc[-1, 1]
        self.translation_table.iloc[:,0] = self.translation_table.apply(lambda x: x.iloc[0].replace("*", "(star)"), axis = 1)
        self.accession = self.translation_table.iloc[-1, 3]
        self.max = self.translation_table.iloc[:,5].dropna().max() + int(self.config.VARIANT_QUERY_PARAMETERS["5p_offset"])
        self.min = self.translation_table.iloc[:,4].dropna().min() - int(self.config.VARIANT_QUERY_PARAMETERS["3p_offset"])
        self.chromosome = self.config.CHROMOSOME_ACCESSIONS[self.accession]
        self.variants = vcf.get_range(self.chromosome, self.min, self.max)
        self.translation_table["ID"] = self.translation_table.apply(
            lambda x: f"c{self.chromosome}_{x.iloc[4]}", 
            axis = 1)

    def __str__(self):
        return self.gene

    def __repr__(self):
        return self.gene

                
