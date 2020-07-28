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

import sys

from .haplotype import Haplotype, NoVariantsException
from .gene import AbstractGene
from . import LOGGING

class Subject:

    def __init__(self, prefix: str, genes: [AbstractGene], config = None) -> None:
        """
        Subject object - manages data and functions for a single sample in a VCF file
        
        Args:
            prefix (str): Subject ID (comes from the VCF file)
            genes ([Gene]): List of gene.Gene objects
        """
        self.config = config
        self.prefix = prefix
        self.called_haplotypes = {}
        for gene in genes:
            try:
                haplotype = Haplotype(gene, self.prefix, self.config)
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
    


