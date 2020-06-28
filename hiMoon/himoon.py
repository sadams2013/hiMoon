from .vcf import VarFile
from .gene import AbstractGene
from .haplotype import Haplotype
from .subject import Subject
from . import CONFIG, set_config

def get_haps_from_vcf(translation_table_path, vcf_file_path, sample_id, config_path = None):
    if config_path:
        set_config(config_path)
    vcf = VarFile(vcf_file_path, sample_id)
    gene = AbstractGene(translation_table_path, vcf = vcf)
    haplotype = Haplotype(gene, sample_id)
    haplotype.table_matcher()
    return haplotype.optimize_hap()

def get_haps_from_variants(translation_table_path, vcf_data, sample_id, config_path = None):
    if config_path:
        set_config(config_path)
    gene = AbstractGene(translation_table_path, variants = vcf_data)
    haplotype = Haplotype(gene, sample_id)
    haplotype.table_matcher()
    return haplotype.optimize_hap()
