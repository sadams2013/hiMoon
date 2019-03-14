"""
Allele picking algorithm.
For a given gene, iterates over haplotypes, finds 100%
matches, and finds the pair that maximizes the number
of matched variations. 
"""

import copy
from itertools import groupby

from . import logging


class AllelePicker:
    def __init__(self, gene: object, haplotypes: list):
        self.haplotypes = haplotypes
        self.gene = gene

    def pick_stars(self) -> tuple:
        """
        Iterates over haplotypes, chooses the most likely pair(s)
        Also returns T/F of if a novel haplotype is likely
        (Or a haplotype not included in the argued translation table)
        """
        hap_sum = lambda x: sum(
            [item for key, item in x.items()]
        )  # Count number of matching haploid alleles
        possible_novel = (
            False
        )  # Set flag for if the subject has a probable novel allele
        ref = Ref(
            self.gene.ref
        )  # Create ref allele object with functions matching haplotype objects
        haps = []  # Init list of haps that match
        var_alleles = []
        for hapname, hapdata in self.haplotypes.items():
            var_alleles += (
                hapdata.positions_matching.keys()
            )  # Add all positions that matched, even if not 100% identity
            if hapdata.match():  # 100% match gets appended to possible haplotypes
                haps.append(hapdata)
        if len(haps) == 0:  # Only ref alleles because no alt alleles matched 100%
            dip_set = [[ref, ref]]
            called_diplotype = [ref, ref]
            if len(var_alleles) > 0:
                possible_novel = True
        elif len(haps) == 1:  # Must be heterozygote or homozygote for a single allele
            called_diplotype = [haps[0], ref]  # Default call is var/ref
            if hap_sum(haps[0].positions_matching) == 2 * len(
                haps[0].positions_matching.keys()
            ):
                # Are there 2 of each allele? change to homozygote var
                called_diplotype = [haps[0], haps[0]]
            elif hap_sum(haps[0].positions_matching) > len(
                haps[0].positions_matching.keys()
            ):
                # More than heterozygote at some markers, but not all of them
                possible_novel = True
            dip_set = [called_diplotype]
            if len(set(var_alleles)) > len(haps[0].positions_matching.keys()):
                # Are there unaccounted for variants in other haplotypes?
                possible_novel = True
        else:  # Must have two non-ref alleles
            possibles = [
                hap for hap in list(set(haps))
            ]  # Eliminate duplicates TODO: Is this needed or just a sanity check?
            possibles2 = copy.copy(
                possibles
            )  # Make copy (needed because we modify in place)
            possible_dips = []  # Init list of possible diplotypes
            for possible in possibles:  # Iterate over all possible individual alleles
                possible_dips += self._pick_best_match(possible, possibles2, ref)
            markers_matched_int = [
                list(dip[0].positions_matching.keys())
                + list(dip[1].positions_matching.keys())
                for dip in possible_dips
            ]  # What markers are matched?
            markers_matched = set(
                [marker for sublist in markers_matched_int for marker in sublist]
            )  # Remove duplicate markers
            max_score = max(
                [
                    len(dip[0].positions_matching) + len(dip[1].positions_matching)
                    for dip in possible_dips
                ]
            )  # Find the maximum possible score in all matched dips
            possible_dips_final = [
                dip
                for dip in possible_dips
                if (len(dip[0].positions_matching) + len(dip[1].positions_matching))
                == max_score
            ]  # Match all dips that meet the maximum score
            dip_set = self.dip_set(possible_dips_final)
            called_diplotype = dip_set[0]
            if len(set(var_alleles)) > max_score:
                # Are there unaccounted for variant alleles?
                possible_novel = True
            for dip in dip_set:
                # Iterate through all possible dips that match
                # Find out if any of the non-shared variants
                # are homozygote, signaling that there are
                # unused variant alleles.
                dl1, dl2 = self._find_unique(
                    dip[0].positions_matching.keys(), dip[1].positions_matching.keys()
                )
                for var in dl1:
                    if dip[0].positions_matching[var] > 1:
                        possible_novel = True
                for var in dl2:
                    if dip[1].positions_matching[var] > 1:
                        possible_novel = True
        return (self.gene.version, called_diplotype, dip_set, possible_novel)

    def _pick_best_match(self, chosen: object, haps: list, ref: object) -> list:
        """
        Given a single haplotype and a list of other possible haplotypes
        Return valid diplotypes
        """
        new_haps = []
        for haplotype in haps:
            haplotype_copy = copy.copy(haplotype)
            haplotype_copy.update_match(chosen)
            new_haps.append(haplotype_copy)
            if len([hapdata for hapdata in new_haps if hapdata.match()]) == 0:
                return [[chosen, ref]]
        return [[chosen, hapdata] for hapdata in new_haps if hapdata.match()]

    def _find_unique(self, l1: list, l2: list) -> list:
        """
        Take two lists, return two lists of unique values from the other
        """
        _l1 = [i for i in l1 if i not in l2]
        _l2 = [i for i in l2 if i not in l1]
        return _l1, _l2

    def dip_set(self, dip_list: list) -> list:
        """
        Removes duplicate diplotypes from a list of diplotypes
        """
        dip_dict = {}
        dip_set_list = []
        for dip in dip_list:
            for hap in dip:
                dip_dict[str(hap)] = hap
        dip_str_list = [tuple(sorted([str(hap) for hap in dip])) for dip in dip_list]
        dip_str_nodups = [k for k, v in groupby(sorted(dip_str_list))]
        for dip in dip_str_nodups:
            dip_set = []
            for hap in dip:
                dip_set.append(dip_dict[str(hap)])
            dip_set_list.append(dip_set)
        return dip_set_list


class Ref:
    """
    Defines the reference allele as an object with dummy methods
    that match the MatchCounter object.
    """

    def __init__(self, name: str):
        self.name = name
        self.positions_matching = {}
        self.rsIDs_matched = []

    def __repr__(self):
        return self.name

    def __str__(self):
        return self.name

