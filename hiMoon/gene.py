import pandas as pd
import sys

from . import logging
from .vcf import VarFile
from .config import ConfigData

class Gene:

    def __init__(self, translation_table: str, config: ConfigData, vcf: VarFile = None, variants = None) -> None:
        """Create a Gene object
        
        Args:
            translation_table (str): path to translation table
            config (Config): config object from config.Config
            vcf (VarFile): parsed VCF object from vcf.VarFile
        """
        self.gene = None
        self.accession = None
        self.config = config
        self.version = self.get_version(translation_table)
        self.translation_table, self.chromosome, self.reference = self.read_translation_table(translation_table, config)
        self.gene = self.translation_table.iloc[-1, 1]
        self.max = self.translation_table.iloc[:,5].dropna().max() + int(config.VARIANT_QUERY_PARAMETERS["5p_offset"])
        self.min = self.translation_table.iloc[:,4].dropna().min() - int(config.VARIANT_QUERY_PARAMETERS["3p_offset"])
        if vcf:
            self.variants = vcf.get_range(self.chromosome, self.min, self.max)
        else:
            self.variants = variants

    def __str__(self):
        return self.gene

    def __repr__(self):
        return self.gene
    
    @staticmethod
    def get_version(translation_table: str) -> str:
        """Get the translation table version
        
        Args:
            translation_table (str): path to translation table file
        
        Returns:
            str: version string from top line of file
        """
        with open(translation_table, 'rt') as trans_file:
            version = trans_file.readline().strip("#version=\n\t")
        return(version)
    
    @staticmethod
    def read_translation_table(translation_table: str, config: ConfigData) -> pd.DataFrame:
        """Read and process a translation table
        
        Args:
            translation_table (str): path to translation table file
            config (ConfigData): config.Config object
        
        Returns:
            pd.DataFrame: translation table as pandas data frame
        """
        translation_table = pd.read_csv(
                                        translation_table, 
                                        skipinitialspace = True,
                                        delim_whitespace=True,
                                        skiprows = 2,
                                        na_values= {4: ".", 5: "."},
                                        dtype = {4: pd.Int64Dtype(), 5: pd.Int64Dtype()},
                                        names = ["Haplotype Name", "Gene", 
                                                "rsID", "ReferenceSequence",
                                                "Variant Start", "Variant Stop",
                                                "Reference Allele", "Variant Allele",
                                                "Type"]
                                        )
        translation_table.iloc[:,0] = translation_table.apply(lambda x: x.iloc[0].replace("*", "(star)"), axis = 1)
        accession = translation_table.iloc[-1, 3]
        chromosome = config.CHROMOSOME_ACCESSIONS[accession]
        translation_table["ID"] = translation_table.apply(lambda x: f"c{chromosome}_{x.iloc[4]}", axis = 1)
        try:
            reference = translation_table[translation_table["ReferenceSequence"] == "REFERENCE"]["Haplotype Name"][0]
        except IndexError:
            reference = "REF"
        return translation_table, chromosome, reference





                
