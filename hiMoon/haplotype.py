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

import sys
import itertools

import pandas as pd
import numpy as np

from pulp import *
from .gene import AbstractGene
from . import LOGGING

class NoVariantsException(Exception):
    """
    Exception to call if a sample is attempted that has zero variants defined. 
    """
    pass

class Haplotype:

    def __init__(self, gene: AbstractGene, sample_prefix: str, config = None) -> None:
        """
        Create a new haplotype object
        This object is not a subclass, but inherits data from the Gene class
        Conceptually, Gene is a fairly abstract class that has meta information used
        by the subject and haplotype classes. 

        Args:
            gene (Gene): gene.Gene object
            sample_prefix (str): Sample ID 
        """
        self.phased = gene.phased
        self.config = config
        self.solver = gene.solver
        self.matched = False
        self.sample_prefix = sample_prefix
        self.genotypes = gene.get_sample_vars(sample_prefix)
        if len(self.genotypes) == 0:
            raise NoVariantsException
        self.translation_table = gene.get_translation_table_copy()
        self.chromosome = gene.chromosome
        self.version = gene.version
        self.reference = gene.reference
    
    def table_matcher(self) -> None:
        """
        Matches variants in the translation table with the subject's variants
        """
        self.matched = True
        matches = self.translation_table.apply(self._match, axis = 1)
        self.translation_table["MATCH"] = [m[0] for m in matches]
        self.translation_table["STRAND"] = [m[1] for m in matches]
        self.translation_table["VAR_ID"] = self.translation_table.apply(
                lambda x: f'{x["ID"]}_{str(x.iloc[7]).strip("<>")}',
                axis = 1
                )
        self.translation_table = self.translation_table.drop(self.translation_table.index[self.translation_table["MATCH"] == 99].tolist())
        haps = self.translation_table["Haplotype Name"]
        no_match = self.translation_table[self.translation_table["MATCH"] == 0].iloc[:,0] # Haplotypes where there is any variant not matching
        drops = []
        for i in no_match.unique():
            if sum([i == k for k in no_match]) > 0:
                drops.append(i)
        self.translation_table = self.translation_table[~self.translation_table.iloc[:,0].isin(drops)] # Drop haplotypes that don't match 100%
        self.variants = self.translation_table.loc[:,["VAR_ID", "MATCH", "STRAND", "Type", "Variant Start"]].drop_duplicates() # List of matched variants
        self.haplotypes = [hap for hap in self.translation_table.iloc[:,0].unique().tolist()] # List of possible haplotypes

    def _mod_vcf_record(self, alt: str, ref: str) -> str:
        """
        Modifies record from VCF to standardized form
        
        Args:
            alt (str): alt allele
            ref (str): ref allele
        
        Returns:
            str: reformatted alt allele
        """
        if alt is None:
            return "-"
        if "<" in alt:
            return f"s{alt.strip('<>')}"
        elif len(ref) > len(alt):
            return "id-"
        elif len(ref) < len(alt):
            return f'id{alt[1:]}' # Remove first position
        else:
            return f's{alt}'

    def _mod_tt_record(self, var_type: str, alt: str) -> list:
        """
        Modifies the translation table ref to a standardized form
        
        Args:
            var_type (str): insertion, deletion, or substitution
            alt (str): allele from translation table
        
        Returns:
            [list]: modified allele as list based on iupac
        """
        alt = alt.strip("<>")
        if var_type == "insertion":
            return [f'id{alt}']
        elif var_type == "deletion":
            return [f'id-']
        else:
            try:
                return [f's{a}' for a in self.config.IUPAC_CODES[alt]]
            except KeyError:
                return [f's{alt}']

    def _match(self, row: pd.core.series.Series) -> (int, int):
        """
        Evaluate match in a single translation table row with a sample

        Args:
            row (pd.core.series.Series): single row from translation table
            genotypes ([type]): list of genotypes

        Returns:
            int: 99 (missing), 0, 1, or 2 (corresponds to the number of matched alleles for a particular position)
        """
        strand = 0
        if row.iloc[8] in ["insertion", "deletion"]:
            new_pos = int(row["ID"].split("_")[1]) - 1
            ID = f'{row["ID"].split("_")[0]}_{new_pos}_SID'
        else:
            ID = row["ID"]
        try:
            genotype = self.genotypes[ID]
        except KeyError: # Not in VCF
            return int(self.config.MISSING_DATA_PARAMETERS["missing_variants"]), strand
        try:
            vcf_geno = [self._mod_vcf_record(g, genotype["ref"]) for g in genotype["alleles"]]
        except AttributeError:
            return int(self.config.MISSING_DATA_PARAMETERS["missing_variants"]), strand
        if vcf_geno == ["-", "-"]:
            return int(self.config.MISSING_DATA_PARAMETERS["missing_variants"]), strand
        tt_alt_geno = self._mod_tt_record(row.iloc[8], row.iloc[7])
        alt_matches = sum([vcf_geno.count(a) for a in tt_alt_geno])
        if alt_matches == 1 and genotype["phased"]:
            strand = 1 if max([vcf_geno.index(a) for a in tt_alt_geno]) == 1 else -1
        elif alt_matches == 2 and genotype["phased"]:
            strand = 3
        return alt_matches, strand
    
    def _haps_from_prob(self, lp_problem: object) -> tuple:
        """
        Take a optimally solved lp problem
        Produce called haplotypes

        Args:
            lp_problem (object): solved lp problem

        Returns:
            tuple: called haplotypes and associated information
        """
        is_ref = False
        haps = []
        variants = []
        for v in lp_problem.variables():
            if v.varValue:
                if v.varValue > 0:
                    if v.name.split("_")[0] == f'c{self.chromosome}':
                        variants.append(v.name)
                    else:
                        haps.append((v.name, v.varValue))
        if len(haps) == 0:
            called = [self.reference, self.reference]
            is_ref = True
        elif len(haps) == 2:
            called = [haps[0][0], haps[1][0]]
        else:
            called = np.array([np.repeat(i[0], i[1]) for i in haps]).flatten().tolist()
            if len(called) == 1:
                called.append(self.reference)
        return called, variants, len(haps), is_ref
    
    def _solve(self, hap_prob: object) -> object:
        if self.solver == "GLPK":
            hap_prob.solve(GLPK(msg=0))
        else:
            hap_prob.solve(PULP_CBC_CMD(msg=0))

    
    def lp_hap(self) -> tuple:
        """
        Build and run the LP problem

        Returns:
            tuple: list of possible haplotypes and list of associated variants
        """
        possible_haplotypes = []
        haplotype_variants = []
        num_vars = self.variants.shape[0]
        num_haps = len(self.haplotypes)
        hap_vars = []
        for hap in self.haplotypes:
            trans = self.translation_table[self.translation_table.iloc[:,0] == hap]
            hap_vars.append([1 if var in trans["VAR_ID"].unique() else 0 for var in self.variants["VAR_ID"]])

        hap_prob = LpProblem("Haplotype Optimization", LpMaximize)

        # Define the haplotypes and variants variables
        haplotypes = [LpVariable(hap, cat = "Integer", lowBound=0, upBound=2) for hap in self.haplotypes]
        variants = [LpVariable(var, cat = "Binary") for var in self.variants["VAR_ID"]]
        # Set constraint of two haplotypes selected
        hap_prob += (lpSum(haplotypes[i] for i in range(num_haps)) <= int(self.config.LP_PARAMS["max_haps"])) # Cannot choose more than x haplotypes (may be increased to find novel sub-alleles)
       # Limit alleles that can be chosen based on zygosity
        for i in range(num_vars): # Iterate over every variant
            # A variant allele can only be used once per haplotype, up to two alleles per variant
            hap_prob += (variants[i] <= (lpSum(hap_vars[k][i] * haplotypes[k] for k in range(num_haps))))
            # A given variant cannot be used more than "MATCH"
            hap_prob += ((lpSum(hap_vars[k][i] * haplotypes[k] for k in range(num_haps))) <= self.variants.iloc[i,1] * variants[i])
            # Any CNV variants defined, if matched with a haplotype, MUST be used
            # Otherwise, variants like CYP2D6*5 will be missed by the other methods
            if self.variants.iloc[i,3] == "CNV":
                hap_prob += ((lpSum(hap_vars[k][i] * haplotypes[k] for k in range(num_haps))) == self.variants.iloc[i,1])
        if self.phased:
            for i in range(num_haps):
                hap_prob += lpSum(haplotypes[i] * self._get_strand_constraint(i, []).size) <= 1 # max one strand
            hap_prob += lpSum(haplotypes[i] * self._get_strand_constraint(i, [0])[0] for i in range(num_haps)) <= 1
            hap_prob += lpSum(haplotypes[i] * self._get_strand_constraint(i, [0])[0] for i in range(num_haps)) >= -1
        # Set to maximize the number of variant alleles used
        hap_prob += lpSum(
            self.translation_table[
                    (self.translation_table.iloc[:,0] == self.haplotypes[i]) &
                    (self.translation_table["MATCH"] > 0)
                ].shape[0] * haplotypes[i] for i in range(num_haps))
        self._solve(hap_prob)
        if hap_prob.status != 1:
            if self.phased:
                LOGGING.warning(f"No feasible solution found, {self.sample_prefix} will be re-attempted with phasing off.")
                self.phased = False
                return None, None
            else:
                LOGGING.warning(f"No feasible solution found, {self.sample_prefix} will not be called")
                return [], []
        else:
            called, variants, hap_len, is_ref = self._haps_from_prob(hap_prob)
            if is_ref:
                possible_haplotypes.append(tuple(called))
                haplotype_variants.append(tuple(variants))
                return possible_haplotypes, haplotype_variants
            max_opt = hap_prob.objective.value()
            opt = max_opt
            while opt >= (max_opt - float(self.config.LP_PARAMS["optimal_decay"])) and not is_ref and hap_prob.status >= 0:
                possible_haplotypes.append(tuple(sorted(called)))
                haplotype_variants.append(tuple(sorted(variants)))
                hap_prob += lpSum([h.value() * h for h in haplotypes]) <= hap_len - 1
                self._solve(hap_prob)
                if hap_prob.status != 1:
                    break
                opt = hap_prob.objective.value()
                new_called, variants, hap_len, is_ref = self._haps_from_prob(hap_prob)
                if new_called == called or len(new_called) == 0:
                    break
                called = new_called

            return possible_haplotypes, haplotype_variants
    

    def _get_strand_constraint(self, i: int, default: list) -> list:
        """
        Helps to assemble the constraint for phased data
        Looks at all strands that are part of a haplotype
        Removes homozygous calls
        Args:
            i (int): haplotype index
            default (list): default return if nothing matches or if all are homozygous

        Returns:
            list: [description]
        """
        sc = self.translation_table[self.translation_table.iloc[:,0] == self.haplotypes[i]]["STRAND"].unique()
        sc = np.delete(sc, np.where(sc == [3]))
        return(sc if sc.size > 0 else np.array(default))


    def optimize_hap(self) -> ():
        """
        Solve for the most likely diplotype

        Returns:
            (): Results
        """
        if not self.matched:
            print("You need to run the table_matcher function with genotyped before you can optimize")
            sys.exit(1)
        called, variants = self.lp_hap()
        if called is None:
            # Happens when a phased call attempt fails
            called, variants = self.lp_hap()
        if len(called) > 1:
            LOGGING.warning(f"Multiple genotypes possible for {self.sample_prefix}.")
        return called, variants
      
