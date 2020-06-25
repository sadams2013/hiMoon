import unittest
import csv

import pandas as pd

from hiMoon import gene, vcf, subject, config, himoon

VCF = vcf.VarFile("hiMoon/tests/test_files/NA12878_chr22.bcf")
CYP2D6_TABLE = "hiMoon/tests/test_files/CYP2D6.NC_000022.11.haplotypes.tsv"
GENE = gene.Gene(CYP2D6_TABLE, config = config.ConfigData(), vcf = VCF)
SUBJ = subject.Subject("NA12878", genes = [GENE])

class TestConfig(unittest.TestCase):

    def test_load(self):
        conf = config.ConfigData("hiMoon/tests/test_files/config.ini")
        assert conf.CHROMOSOME_ACCESSIONS["TESTCHROM"] == "testvalue"
    
    def test_dummy_file(self):
        conf = config.ConfigData()
        assert conf.CHROMOSOME_ACCESSIONS["NC_000001.11"] == "1"

class TestGene(unittest.TestCase):

    def test_gene(self):
        assert GENE.gene == "CYP2D6"
    
    def test_translation_table_version(self):
        assert gene.Gene.get_version(CYP2D6_TABLE) == "pharmvar-4.1.5"
    
    def test_translation_table_read(self):
        table = gene.Gene.read_translation_table(CYP2D6_TABLE, config.ConfigData())[0]
        assert type(table) == pd.DataFrame


class TestSubject(unittest.TestCase):

    def test_subject_prefix(self):
        self.assertEqual(SUBJ.prefix, "NA12878")
    
    def test_called_haplotypes(self):
        HAPS = sorted([i.split(".")[0] for i in SUBJ.called_haplotypes["CYP2D6"]["HAPS"][1]])
        self.assertEqual(HAPS, sorted(["CYP2D6(star)3", "CYP2D6(star)4"]))

class TestVCF(unittest.TestCase):

    def test_samples(self):
        self.assertEqual(VCF.samples[0], "NA12878")

class TestHiMoon(unittest.TestCase):

    def test_himoon(self):
        haps = himoon.get_haps_from_vcf(
            "hiMoon/tests/test_files/CYP2D6.NC_000022.11.haplotypes.tsv",
            "hiMoon/tests/test_files/NA12878_chr22.bcf",
            "NA12878"
        )
        HAPS = sorted([i.split(".")[0] for i in haps[1]])
        self.assertEqual(HAPS, sorted(["CYP2D6(star)3", "CYP2D6(star)4"]))
