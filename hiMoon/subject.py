"""
Defines the Subject class
This serves to coordinate multiple gene lookups
for the same subject.
"""

import sys

from .haplotype import Haplotype
from .gene import Gene
from . import logging

class Subject:
    """
    Describes a subject object
    Subject have called haplotypes and haplotype match objects
    These contain verbose information about how a haplotype was called
    """

    def __init__(self, prefix: str, genes: list, config):
        self.prefix = prefix
        self.called_haplotypes = {}
        self.config = config
        for gene in genes:
            genotypes = self._get_vars(gene)
            haplotype = Haplotype(gene)
            haplotype.table_matcher(genotypes)
            self.called_haplotypes[str(gene)] = haplotype.optimize_hap()
    
    def _get_vars(self, gene):
        genotypes = {}
        for var_id, sub_vars in gene.variants.items():
            try:
                genotypes[var_id] = sub_vars[self.prefix]
            except KeyError:
                pass
        return genotypes

