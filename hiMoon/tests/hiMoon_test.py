import unittest
import csv

from hiMoon import gene, vcf, subject, alignment, cnv


class TestGene(unittest.TestCase):

    def test_gene(self):
        self.gene_obj = gene.Gene("hiMoon/tests/CYP2D6.NC_000022.11haplotypes.tsv")
        self.assertEqual(self.gene_obj.gene, "CYP2D6", "Make sure that gene name is correct")
        self.assertEqual(self.gene_obj.position_max(), 42133392, "Ensure max value is correctly called")
        self.assertEqual(self.gene_obj.position_min(), 42125310, "Ensure min value is correctly called")
        self.assertEqual(self.gene_obj.chromosome(), 22, "Ensure chromosome is correctly called")
        self.assertEqual(self.gene_obj.version, "pharmvar-3.3", "Ensure version is correctly called")

class TestVCF(unittest.TestCase):

    def test_vcf(self):
        self.vcf_obj = vcf.VCFParse("hiMoon/tests/test.cyp2d6.vcf.gz")


class TestSubject(unittest.TestCase):

    def test_subject(self):
        self.subject_object_vcf = subject.Subject(
            prefix = "",
            genes = [gene.Gene("hiMoon/tests/CYP2D6.NC_000022.11haplotypes.tsv")],
            vcf_file = "hiMoon/tests/test.cyp2d6.vcf.gz"
        )
        self.subject_object_bam = subject.Subject(
            prefix = "",
            genes = [gene.Gene("hiMoon/tests/CYP2D6.NC_000022.11haplotypes.tsv")],
            alignment_file = "hiMoon/tests/NA18565.cyp2d6.bam"
        )

class TestAlignment(unittest.TestCase):

    def test_alignment(self):
        self.alignment_obj = alignment.AlignmentData(
            alignment_file = "hiMoon/tests/NA18565.cyp2d6.bam"
            )
        self.var_range_1 = self.alignment_obj.get_range(22, 42128240, 42128245)
        self.var_range_2 = self.alignment_obj.get_range(22, 42110000, 42140000)
    
    def test_alignment_insertion(self):
        alignment_obj = alignment.AlignmentData(
            alignment_file = "hiMoon/tests/NA19239.cyp2d6.bam"
            )
        alleles = alignment_obj.call_alleles("22", 42130650, 42130660)

class TestCNV(unittest.TestCase):

    def test_cnv(self):
          case_bam = alignment.AlignmentData(
               alignment_file="hiMoon/tests/NA18945.cyp2d6.bam"
          )
          control_bam = alignment.AlignmentData(
               alignment_file="hiMoon/tests/NA19239.cyp2d6.bam"
          )
          case = case_bam.get_read_array("22", 42128000,42130000)
          control = control_bam.get_read_array("22", 42128000,42130000)
          a = cnv.find_cnv(case, control, 1000)
