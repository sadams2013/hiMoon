import unittest
import csv

from hiMoon import gene, vcf, subject, config

config = config.ConfigData("hiMoon/tests/default_config.ini")

class TestGene(unittest.TestCase):

    def test_gene(self):
        self.gene_obj = gene.Gene("hiMoon/tests/CYP2D6.NC_000022.11haplotypes.tsv", config = config)


class TestSubjet(unittest.TestCase):

    def test_subject(self):
        self.gene_obj = gene.Gene("hiMoon/tests/CYP2D6.NC_000022.11haplotypes.tsv", config = config)
        subj = subject.Subject(
            "test",
            genes = [self.gene_obj],
            config = config
            )

