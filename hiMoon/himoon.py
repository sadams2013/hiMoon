from .vcf import VarFile
from .gene import Gene
from .config import ConfigData
from .haplotype import Haplotype
from .subject import Subject

def get_haps(translation_table_path, vcf_file_path, sample_id, config_path):
    vcf = VarFile(vcf_file_path, sample_id)
    config = ConfigData(config_path)
    gene = Gene(translation_table_path, config, vcf)
    haplotype = Haplotype(gene, sample_id)
    subject = Subject(sample_id, genes = [gene])
    return subject.called_haplotypes, gene.translation_table, gene.variants[sample_id]

