import pandas as pd
import numpy as np
import sys

from pulp import *


class Haplotype:

    def __init__(self, gene):
        self.translation_table = gene.translation_table
        self.ref = gene.ref
        self.chromosome = gene.chromosome
        self.version = gene.version

    def table_matcher(self, genotypes):
        self.matched = True
        self.translation_table["MATCH"] = self.translation_table.apply(
            lambda x: self._match(x, genotypes),
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

    def _mod_alt(self, alt, ref):
        # if its a del, needs to return -s
        # if its an ins, needs to return just what is inserted
        if len(ref) > len(alt):
            return "id-"
        elif len(ref) > 1:
            return f'id{alt[1:]}' # Remove first position
        else:
            return f's{alt}'
    
    def _mod_ref(self, var_type, alt):
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
        geno = [self._mod_alt(g, genotype["ref"]) for g in genotype["alleles"]]
        tt_alt = self._mod_ref(row.iloc[8], row.iloc[7])
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
                        variants.append(f"{v.name}_{v.varValue}")
                    else:
                        haps.append((v.name, v.varValue))
        if len(haps) == 0:
            called = ["ref", "ref"]
        else:
            called = np.array([np.repeat(i[0], i[1]) for i in haps]).flatten().tolist()
            if len(called) == 1:
                called.append("ref")
        return(self.version, haps, variants, called)