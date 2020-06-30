from .vcf import VarFile
from .gene import AbstractGene
from .haplotype import Haplotype
from .subject import Subject
from . import CONFIG, set_config

def get_haps_from_vcf(translation_table_path: str, vcf_file_path: str, sample_id: str, config_path = None) -> tuple:
    """
    Provide a VCF file, sample ID, and translation table
    Get called haplotypes and additional information

    Args:
        translation_table_path (str): [description]
        vcf_file_path (str): [description]
        sample_id (str): [description]
        config_path ([type], optional): [description]. Defaults to None.

    Returns:
        tuple: translation_table_version, called_haplotypes, variants_associated_with_haplotye, matched_translation_table
    """
    if config_path:
        set_config(config_path)
    vcf = VarFile(vcf_file_path, sample_id)
    gene = AbstractGene(translation_table_path, vcf = vcf)
    haplotype = Haplotype(gene, sample_id)
    haplotype.table_matcher()
    return haplotype.optimize_hap()

def get_haps_from_variants(translation_table_path: str, vcf_data: str, sample_id: str, config_path = None) -> tuple:
    """
    Same as get_haps_from_vcf, but bypasses the VCF file so that you can provide formatted variants from another input
    Get called haplotypes and additional information

    Args:
        translation_table_path (str): [description]
        vcf_file_path (str): [description]
        sample_id (str): [description]
        config_path ([type], optional): [description]. Defaults to None.

    Returns:
        tuple: translation_table_version, called_haplotypes, variants_associated_with_haplotye, matched_translation_table
    """
    if config_path:
        set_config(config_path)
    gene = AbstractGene(translation_table_path, variants = vcf_data)
    haplotype = Haplotype(gene, sample_id)
    haplotype.table_matcher()
    return haplotype.optimize_hap()
