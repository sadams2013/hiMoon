import sys

from .haplotype import Haplotype
from .gene import Gene
from . import logging

class Subject:

    def __init__(self, prefix: str, genes: [Gene]) -> None:
        """Subject object - manages data and functions for a single sample in a VCF file
        
        Args:
            prefix (str): Subject ID (comes from the VCF file)
            genes ([Gene]): List of gene.Gene objects
        """
        self.prefix = prefix
        self.called_haplotypes = {}
        for gene in genes:
            haplotype = Haplotype(gene, self.prefix)
            haplotype.table_matcher()
            self.called_haplotypes[str(gene)] = haplotype.optimize_hap()
    
    def __str__(self):
        return self.prefix

    def __repr__(self):
        return self.prefix
    


