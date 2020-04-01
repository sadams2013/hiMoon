import sys

from .haplotype import Haplotype
from .gene import Gene
from .config import ConfigData
from . import logging

class Subject:

    def __init__(self, prefix: str, genes: [Gene], config: ConfigData) -> None:
        """Subject object - manages data and functions for a single sample in a VCF file
        
        Args:
            prefix (str): Subject ID (comes from the VCF file)
            genes ([Gene]): List of gene.Gene objects
            config (ConfigData): config.Config object
        """
        self.prefix = prefix
        self.called_haplotypes = {}
        self.config = config
        for gene in genes:
            genotypes = self._get_vars(gene)
            haplotype = Haplotype(gene, genotypes)
            haplotype.table_matcher()
            self.called_haplotypes[str(gene)] = haplotype.optimize_hap()
    
    def __str__(self):
        return self.prefix

    def __repr__(self):
        return self.prefix
    
    def _get_vars(self, gene: Gene) -> dict:
        """The Gene object has the parsed VCF with all samples. This gets the variants for just this subject.
        
        Args:
            gene (Gene): Gene object
        
        Returns:
            dict: dictionary of variants for this subject
        """
        genotypes = {}
        for var_id, sub_vars in gene.variants.items():
            try:
                genotypes[var_id] = sub_vars[self.prefix]
            except KeyError:
                pass
        return genotypes

