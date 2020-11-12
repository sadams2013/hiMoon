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

import pandas as pd
import sys

from . import LOGGING
from .vcf import VarFile

class AbstractGene:
    """
    Abstract gene class, conains top level information that is available to muliple
    names haplotypes and subjects. 
    """

    def __init__(self, translation_table: str, vcf: VarFile = None, 
                    variants = None, solver: str = "CBC",
                    config = None, phased = False) -> None:
        """
        Create a Gene object
        
        Args:
            translation_table (str): path to translation table
            vcf (VarFile): parsed VCF object from vcf.VarFile
        """
        self.config = config
        self.phased = phased
        self.solver = solver
        self.gene = None
        self.accession = None
        self.read_translation_table(translation_table)
        self.gene = self.translation_table.iloc[-1, 1]
        self.max = self.translation_table.iloc[:,5].dropna().max() + int(self.config.VARIANT_QUERY_PARAMETERS["5p_offset"])
        self.min = self.translation_table.iloc[:,4].dropna().min() - int(self.config.VARIANT_QUERY_PARAMETERS["3p_offset"])
        if vcf:
            self.variants = vcf.get_range(self.chromosome, self.min, self.max)
        else:
            self.variants = variants

    def __str__(self):
        return self.gene

    def __repr__(self):
        return self.gene
    
    def get_sample_vars(self, sample: str) -> dict:
        """
        The gene contains variants for all samples in the VCF
        This function parses variants for a single sample. 

        Args:
            sample (str): sample ID

        Returns:
            dict: single sample variants from VCF
        """
        sample_vars = {}
        for var_id, sub_vars in self.variants.items():
            try:
                sample_vars[var_id] = sub_vars[sample]
            except KeyError:
                pass
        return sample_vars
    
    def get_translation_table_copy(self) -> pd.DataFrame:
        """
        Get a deep copy of a translation table that can be modified
        for a single subject

        Returns:
            pd.DataFrame: deep copy of the associated translation table
        """
        return(self.translation_table.copy(deep = True))
    
    def _merge_tables(self, translation_table, cnv_table):
        new_rows = []
        cnv_table["BASE"] = cnv_table.apply(lambda x: x["Haplotype Name"].split("_")[0], axis = 1)
        cnv_table["SUFFIX"] = cnv_table.apply(lambda x: x["Haplotype Name"].split("_")[-1], axis = 1)
        translation_table["BASE"] = translation_table.apply(lambda x: x["Haplotype Name"].split(".")[0], axis = 1)
        translation_table["SUFFIX"] = translation_table.apply(lambda x: x["Haplotype Name"].split(".")[-1], axis = 1)
        for index, row in cnv_table.iterrows():
            trans_base = translation_table[translation_table["BASE"] == row["BASE"]]
            trans_suffixes = trans_base["SUFFIX"]
            for suf in trans_suffixes:
                new_rows.append(
                    {
                        "Haplotype Name": f'{row["BASE"]}.{suf}_{row["SUFFIX"]}',
                        "Gene": row["Gene"],
                        "rsID": row["rsID"],
                        "ReferenceSequence": row["ReferenceSequence"],
                        "Variant Start": row["Variant Start"],
                        "Variant Stop": row["Variant Stop"],
                        "Reference Allele": row["Reference Allele"],
                        "Variant Allele": row["Variant Allele"],
                        "Type": row["Type"]
                    }
                )
            for i, r in trans_base.iterrows():
                new_rows.append(
                    {
                        "Haplotype Name": f'{row["BASE"]}.{r["SUFFIX"]}_{row["SUFFIX"]}',
                        "Gene": r["Gene"],
                        "rsID": r["rsID"],
                        "ReferenceSequence": r["ReferenceSequence"],
                        "Variant Start": r["Variant Start"],
                        "Variant Stop": r["Variant Stop"],
                        "Reference Allele": r["Reference Allele"],
                        "Variant Allele": r["Variant Allele"],
                        "Type": r["Type"]
                    }
                )
        new_table = translation_table.drop(["BASE", "SUFFIX"], axis = 1)
        added_rows = pd.DataFrame(new_rows)
        return pd.concat([new_table, added_rows, cnv_table], ignore_index=True)

    def get_type(self, vtype: str) -> str:
        """
        Simple helper function to assign a CNV type if a variant is not SID

        Args:
            vtype (str): SID or CNV

        Returns:
            str: CNV if CNV, else SID
        """
        return "CNV" if vtype == "CNV" else "SID"

    def read_translation_table(self, translation_table: str) -> pd.DataFrame:
        """Read and process a translation table
        
        Args:
            translation_table (str): path to translation table file
            config (ConfigData): config.Config object
        
        Returns:
            pd.DataFrame: translation table as pandas data frame
        """
        with open(translation_table, 'rt') as trans_file:
            self.version = trans_file.readline().strip("#version=\n\t")
        self.translation_table = pd.read_csv(
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
        try:
            cnv_table = pd.read_csv(
                                    translation_table.replace(".tsv", ".cnv"),
                                    skipinitialspace = True,
                                    delim_whitespace=True,
                                    skiprows = 2,
                                    na_values= {4: ".", 5: "."},
                                    dtype = {4: pd.Int64Dtype(), 5: pd.Int64Dtype(), 6: pd.StringDtype(), 7: pd.StringDtype()},
                                    names = ["Haplotype Name", "Gene", 
                                            "rsID", "ReferenceSequence",
                                            "Variant Start", "Variant Stop",
                                            "Reference Allele", "Variant Allele",
                                            "Type"]
            )
        except FileNotFoundError:
            LOGGING.info("No CNV definition file found, proceeding with only SID variants.")
            cnv_table = None
        if cnv_table is not None:
            self.translation_table = self._merge_tables(self.translation_table, cnv_table)
        self.translation_table.iloc[:,0] = self.translation_table.apply(lambda x: x.iloc[0].replace("*", "(star)"), axis = 1)
        try:
            self.reference = self.translation_table[self.translation_table["rsID"] == "REFERENCE"]["Haplotype Name"][0]
        except IndexError:
            self.reference = "REF"
        self.translation_table = self.translation_table[self.translation_table["ReferenceSequence"] != "."]
        self.accession = self.translation_table.iloc[-1, 3]
        self.chromosome = self.config.CHROMOSOME_ACCESSIONS[self.accession]
        self.translation_table["ID"] = self.translation_table.apply(lambda x: f"c{self.chromosome}_{x['Variant Start']}_{self.get_type(x['Type'])}", axis = 1)
        self.translation_table["EXCLUDE"] = 0
        