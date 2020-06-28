import pandas as pd
import sys

from . import LOGGING, CONFIG
from .vcf import VarFile

class AbstractGene:
    """
    Abstract gene class, conains top level information that is available to muliple
    names haplotypes and subjects. 
    """

    def __init__(self, translation_table: str, vcf: VarFile = None, variants = None) -> None:
        """Create a Gene object
        
        Args:
            translation_table (str): path to translation table
            vcf (VarFile): parsed VCF object from vcf.VarFile
        """
        self.gene = None
        self.accession = None
        self.translation_table, self.chromosome, self.reference = self.read_translation_table(translation_table)
        self.gene = self.translation_table.iloc[-1, 1]
        self.max = self.translation_table.iloc[:,5].dropna().max() + int(CONFIG.VARIANT_QUERY_PARAMETERS["5p_offset"])
        self.min = self.translation_table.iloc[:,4].dropna().min() - int(CONFIG.VARIANT_QUERY_PARAMETERS["3p_offset"])
        if vcf:
            self.variants = vcf.get_range(self.chromosome, self.min, self.max)
        else:
            self.variants = variants

    def __str__(self):
        return self.gene

    def __repr__(self):
        return self.gene
    
    def get_sample_vars(self, sample: str) -> dict:
        """[summary]

        Args:
            sample (str): [description]

        Returns:
            dict: [description]
        """
        sample_vars = {}
        for var_id, sub_vars in self.variants.items():
            try:
                sample_vars[var_id] = sub_vars[sample]
            except KeyError:
                pass
        return sample_vars
    
    def read_translation_table(self, translation_table: str) -> pd.DataFrame:
        """
        Read and process a translation table, returns a data frame.
        This is separated from the class instance such that each subject has
        their own translation table instance that gets modified in haplotype. 
        
        Args:
            translation_table (str): path to translation table file
            CONFIG (ConfigData): CONFIG.Config object
        
        Returns:
            pd.DataFrame: translation table as pandas data frame
        """
        with open(translation_table, 'rt') as trans_file:
            self.version = trans_file.readline().strip("#version=\n\t")
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
        chromosome = CONFIG.CHROMOSOME_ACCESSIONS[accession]
        translation_table["ID"] = translation_table.apply(lambda x: f"c{chromosome}_{x.iloc[4]}", axis = 1)
        try:
            reference = translation_table[translation_table["ReferenceSequence"] == "REFERENCE"]["Haplotype Name"][0]
        except IndexError:
            reference = "REF"
        return translation_table, chromosome, reference
    
    def get_translation_table_copy(self) -> pd.DataFrame:
        """
        Get a deep copy of a translation table that can be modified
        for a single subject

        Returns:
            pd.DataFrame: deep copy of the associated translation table
        """
        return(self.translation_table.copy(deep = True))
