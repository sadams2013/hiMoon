import unittest
import yaml
import glob
import os

from hiMoon import gene, vcf, subject, config, himoon

PATH = os.path.dirname(os.path.abspath(__file__))
CONF = config.ConfigData()

def load_definition_file(file_path: str) -> dict:
    """
    Wrapper for yaml.load

    Args:
        file_path (str): path/to/yml/file

    Returns:
        dict: dict representation of yml
    """
    with open(file_path, "r") as sample_allele_data:
        try:
            return(yaml.load(sample_allele_data, Loader = yaml.CLoader))
        except AttributeError: # Happens if libyaml headers are not available
            return(yaml.load(sample_allele_data, Loader = yaml.Loader))

class TestGenes(unittest.TestCase):

    def test_alleles(self):
        definition_files = glob.glob(PATH + "/test_files/*.yml")
        for definition_file in definition_files:
            sample_alleles = load_definition_file(definition_file)
            vcf_file = vcf.VarFile(PATH + "/test_files/" + sample_alleles["VCF"])
            gene_obj = gene.Gene(PATH + "/test_files/" + sample_alleles["TRANSLATION_TABLE"], CONF, vcf = vcf_file)
            for sample in sample_alleles["SAMPLES"]:
                subj = subject.Subject(sample["ID"], genes = [gene_obj])
                haps = sorted([i.split(".")[0] for i in subj.called_haplotypes[str(gene_obj)]["HAPS"][1]])
                self.assertEqual(haps, sorted(sample["ALLELES"]), f"{sample['ID']} @ {str(gene)} should be {'/'.join(sample['ALLELES'])}. Found {'/'.join(haps)}")




# Test specific genotype calls