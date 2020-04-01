import unittest
import csv

import pandas as pd

from hiMoon import gene, vcf, subject, config

CONFIG = config.ConfigData("hiMoon/tests/config.ini")
VCF = vcf.VarFile("hiMoon/tests/NA12878_chr22.bcf")
CYP2D6_TABLE = "hiMoon/tests/CYP2D6.NC_000022.11.haplotypes.tsv"
GENE = gene.Gene(CYP2D6_TABLE, config = CONFIG, vcf = VCF)
SUBJ = subject.Subject("NA12878",genes = [GENE], config = CONFIG)

class TestGene(unittest.TestCase):

    def test_gene(self):
        assert GENE.gene == "CYP2D6"
    
    def test_translation_table_version(self):
        assert gene.Gene.get_version(CYP2D6_TABLE) == "pharmvar-4.1.4"
    
    def test_translation_table_read(self):
        table = gene.Gene.read_translation_table(CYP2D6_TABLE, CONFIG)[0]
        assert type(table) == pd.DataFrame


class TestSubject(unittest.TestCase):

    def test_subject_prefix(self):
        assert SUBJ.prefix == "NA12878"
    
    def test_called_haplotypes(self):
        assert SUBJ.called_haplotypes["CYP2D6"][1] == [("CYP2D6(star)3", 1.0), ("CYP2D6(star)4.001", 1.0)]

class TestVCF(unittest.TestCase):

    def test_samples(self):
        assert VCF.samples[0] == "NA12878"
