import yaml
import glob
import os
from parameterized import parameterized

from hiMoon import gene, vcf, subject, himoon

PATH = os.path.dirname(os.path.abspath(__file__))

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

def prep_samples() -> []:
    gene_samples = []
    definition_files = glob.glob(PATH + "/test_files/*.yaml")
    for definition_file in definition_files:
        sample_alleles = load_definition_file(definition_file)
        vcf_file = vcf.VarFile(PATH + "/test_files/" + sample_alleles["VCF"])
        gene_obj = gene.AbstractGene(PATH + "/test_files/" + sample_alleles["TRANSLATION_TABLE"], vcf = vcf_file)
        gene_samples += [(subject.Subject(sample["ID"], genes = [gene_obj]), gene_obj, sample) for sample in sample_alleles["SAMPLES"]]
    return(gene_samples)

@parameterized.expand(prep_samples())
def test_sample_haplotypes(subj, gene_obj, sample):
    haps = sorted([i.split(".")[0] for i in subj.called_haplotypes[str(gene_obj)]["HAPS"][1]])
    print(f"{sample['ID']} @ {str(gene_obj)} should be {'/'.join(sample['ALLELES'])}. Found {'/'.join(haps)}")
    assert haps == sorted(sample["ALLELES"])


if __name__ == "__main__":
    test_sample_haplotypes()
    