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

import unittest
import csv
import os

import pandas as pd

from hiMoon import gene, vcf, subject, config, himoon, get_config

CONFIG = get_config()

PATH = os.path.dirname(os.path.abspath(__file__))

VCF = vcf.VarFile(PATH + "/test_files/vcf/test_samples.bcf")
CYP2D6_TABLE = PATH + "/test_files/translation_tables/CYP2D6.NC_000022.11.haplotypes.tsv"
GENE = gene.AbstractGene(CYP2D6_TABLE, vcf = VCF, config = CONFIG)
SUBJ = subject.Subject("NA12878", genes = [GENE], config = CONFIG)

class TestConfig(unittest.TestCase):

    def test_load(self):
        conf = get_config(PATH + "/test_files/config.ini")
        assert conf.CHROMOSOME_ACCESSIONS["TESTCHROM"] == "testvalue"
    
    def test_dummy_file(self):
        assert CONFIG.CHROMOSOME_ACCESSIONS["NC_000001.11"] == "1"

class TestGene(unittest.TestCase):

    def test_gene(self):
        assert GENE.gene == "CYP2D6"
    
    def test_translation_table_version(self):
        assert GENE.version == "pharmvar-4.1.6"
    
    def test_reference(self):
        assert GENE.reference == "CYP2D6(star)1"

class TestSubject(unittest.TestCase):

    def test_subject_prefix(self):
        self.assertEqual(SUBJ.prefix, "NA12878")

class TestVCF(unittest.TestCase):

    def test_samples(self):
        self.assertEqual(VCF.samples[0], "HG00111")
