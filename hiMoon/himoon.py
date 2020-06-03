from .vcf import VarFile
from .gene import Gene
from .config import ConfigData
from .haplotype import Haplotype

def get_haps_from_vcf(translation_table_path, vcf_file_path, sample_id, config_path):
    vcf = VarFile(vcf_file_path, sample_id)
    config = ConfigData(config_path)
    gene = Gene(translation_table_path, config, vcf = vcf)
    haplotype = Haplotype(gene, sample_id)
    haplotype.table_matcher()
    return haplotype.optimize_hap()

def get_haps_from_variants(translation_table_path, vcf_data, sample_id, config_path):
    config = ConfigData(config_path)
    gene = Gene(translation_table_path, config, variants = vcf_data)
    haplotype = Haplotype(gene, sample_id)
    haplotype.table_matcher()
    return haplotype.optimize_hap()