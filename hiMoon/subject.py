import sys

from .haplotype import Haplotype, NoVariantsException
from .gene import AbstractGene
from . import LOGGING

class Subject:

    def __init__(self, prefix: str, genes: [AbstractGene]) -> None:
        """Subject object - manages data and functions for a single sample in a VCF file
        
        Args:
            prefix (str): Subject ID (comes from the VCF file)
            genes ([Gene]): List of gene.Gene objects
        """
        self.prefix = prefix
        self.called_haplotypes = {}
        for gene in genes:
            try:
                haplotype = Haplotype(gene, self.prefix)
                haplotype.table_matcher()
                self.called_haplotypes[str(gene)] = {
                    "HAPS": haplotype.optimize_hap(),
                    "CONTIG": gene.chromosome}
            except NoVariantsException:
                LOGGING.warning(f"{self.prefix} has no variants, returning NA")
                self.called_haplotypes[str(gene)] = {
                    "HAPS": ("NA", "NA", "NA", "NA"),
                    "CONTIG": gene.chromosome}

    
    def __str__(self):
        return self.prefix

    def __repr__(self):
        return self.prefix
    


