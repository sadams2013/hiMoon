import pandas as pd
import numpy as np
import sys

from pulp import *

from .gene import Gene


class Haplotype:

    def __init__(self, gene: Gene, sample_prefix: str) -> None:
        """Subject level haplotypes for a given gene
        
        Args:
            gene (Gene): gene.Gene object
            genotypes (dict): subject's genotypes
        """
        self.translation_table = gene.translation_table
        self.chromosome = gene.chromosome
        self.version = gene.version
        self.sample_prefix = sample_prefix
        self.genotypes = self._get_vars(gene)
        self.reference = gene.reference
        self.iupac = gene.config.IUPAC_CODES
    
    def _get_vars(self, gene) -> dict:
        """The Gene object has the parsed VCF with all samples. This gets the variants for just this subject.
        
        Returns:
            dict: dictionary of variants for this subject
        """
        genotypes = {}
        for var_id, sub_vars in gene.variants.items():
            try:
                genotypes[var_id] = sub_vars[self.sample_prefix]
            except KeyError:
                pass
        return genotypes

    def table_matcher(self) -> None:
        """Matches translation table with the subject's genotypes
        """
        self.matched = True
        self.translation_table["MATCH"] = self.translation_table.apply(
            lambda x: self._match(x, self.genotypes),
            axis = 1
        )
        self.translation_table["VAR_ID"] = self.translation_table.apply(
                lambda x: f'{x["ID"]}_{x.iloc[7]}',
                axis = 1
                )
        self.translation_table = self.translation_table.drop(self.translation_table.index[self.translation_table["MATCH"] == 99].tolist())
        no_match = self.translation_table[self.translation_table["MATCH"] == 0].iloc[:,0].unique() # Haplotypes where there is any variant not matching
        self.translation_table = self.translation_table[~self.translation_table.iloc[:,0].isin(no_match)] # Drop haplotypes that don't match 100%
        self.variants = self.translation_table.loc[:,["VAR_ID", "MATCH"]].drop_duplicates() # List of matched variants
        self.haplotypes = [hap for hap in self.translation_table.iloc[:,0].unique().tolist()] # List of possible haplotypes

    def _mod_vcf_record(self, alt: str, ref: str) -> str:
        """Modifies alt from VCF to standardized form
        
        Args:
            alt (str): alt allele
            ref (str): ref allele
        
        Returns:
            str: reformatted alt allele
        """
        # if its a del, needs to return -s
        # if its an ins, needs to return just what is inserted
        if len(ref) > len(alt):
            return "id-"
        elif len(ref) > 1:
            return f'id{alt[1:]}' # Remove first position
        else:
            return f's{alt}'
    
    def _mod_tt_record(self, var_type: str, alt: str) -> str:
        """Modifies the translation table ref to a standardized form
        
        Args:
            var_type (str): insertion, deletion, or substitution
            alt (str): allele from translation table
        
        Returns:
            [str]: modified allele 
        """
        if var_type == "insertion":
            return f'id{alt}'
        elif var_type == "deletion":
            return f'id-'
        else:
            return f's{alt}'


    def _match(self, row, genotypes):
        if row.iloc[8] in ["insertion", "deletion"]:
            new_pos = int(row["ID"].split("_")[1]) - 1
            ID = f'{row["ID"].split("_")[0]}_{new_pos}'
        else:
            ID = row["ID"]
        try:
            genotype = genotypes[ID]
        except KeyError:
            return 99
        geno = [self._mod_vcf_record(g, genotype["ref"]) for g in genotype["alleles"]]
        tt_alt = self._mod_tt_record(row.iloc[8], row.iloc[7])
        alt_matches = geno.count(tt_alt)
        return(alt_matches)
    
    def optimize_hap(self):
        """
        Goal: Maximize the number of genotypes used
        Constraints: 
            1. Number of haplotypes == 2
            2. ???
        """
        if not self.matched:
            print("You need to run the table_matcher function with genotyped before you can optimize")
            sys.exit(1)

        num_vars = self.variants.shape[0]
        num_haps = len(self.haplotypes)

        hap_vars = []

        for hap in self.haplotypes:
            trans = self.translation_table[self.translation_table.iloc[:,0] == hap]
            hap_vars.append([1 if var in trans["VAR_ID"].unique() else 0 for var in self.variants["VAR_ID"]])

        hap_prob = LpProblem("Haplotype Optimization", LpMaximize)
        
        # Define the haplotypes variable
        haplotypes = [LpVariable(hap, cat = "LpInteger", lowBound=0, upBound=2) for hap in self.haplotypes]
        variants = [LpVariable(var, cat = "Binary") for var in self.variants["VAR_ID"]]
        

        # Set constraint of two haplotypes selected
        hap_prob += (lpSum(haplotypes[i] for i in range(num_haps)) <= 2) # Cannot choose more than two haplotypes

        # Limit alleles that can be chosen based on zygosity
        for i in range(num_vars): # Iterate over every variant
            # A variant can only be used once per haplotype
            hap_prob += (variants[i] <= (lpSum(hap_vars[k][i] * haplotypes[k] for k in range(num_haps))))
            # A given variant cannot be used more than "MATCH"
            hap_prob += ((lpSum(hap_vars[k][i] * haplotypes[k] for k in range(num_haps))) <= self.variants.iloc[i,1] * variants[i])

        # Set to maximize the number of variant alleles used
        hap_prob += lpSum(
            self.translation_table[
                self.translation_table.iloc[:,0] == self.haplotypes[i]
                ]["MATCH"].sum() * haplotypes[i] for i in range(num_haps))

        hap_prob.solve(GLPK(msg=0))
        haps = []
        variants = []

        for v in hap_prob.variables():
            if v.varValue:
                if v.varValue > 0:
                    if v.name.split("_")[0] == f'c{self.chromosome}':
                        variants.append(v.name)
                    else:
                        haps.append((v.name, v.varValue))
        if len(haps) == 0:
            called = [self.reference, self.reference]
        elif len(haps) == 2:
            called = [haps[0][0], haps[1][0]]
        else:
            called = np.array([np.repeat(i[0], i[1]) for i in haps]).flatten().tolist()
            if len(called) == 1:
                called.append(self.reference)
        return self.version, called, variants
        