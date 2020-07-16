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

from .vcf import VarFile
from .gene import AbstractGene
from .haplotype import Haplotype
from .subject import Subject
from . import get_config

def get_haps_from_vcf(translation_table_path: str, vcf_file_path: str, 
                        sample_id: str, solver: str = "CBC", 
                        config_path: str = None, phased = False) -> tuple:
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
    config = get_config(config_path)
    vcf = VarFile(vcf_file_path, sample_id)
    gene = AbstractGene(translation_table_path, vcf = vcf, solver = solver, config = config, phased = phased)
    haplotype = Haplotype(gene, sample_id, config = config)
    haplotype.table_matcher()
    return haplotype.optimize_hap()

def get_haps_from_variants(translation_table_path: str, vcf_data: str, 
                            sample_id: str, solver: str = "CBC", 
                            config_path: str = None, phased = False) -> tuple:
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
    config = get_config(config_path)
    gene = AbstractGene(translation_table_path, variants = vcf_data, solver = solver, config = config, phased = phased)
    haplotype = Haplotype(gene, sample_id, config = config)
    haplotype.table_matcher()
    return haplotype.optimize_hap()
