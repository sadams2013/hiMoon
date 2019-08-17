"""
Defines the MatchCounter class, which is an abstraction of haplotype
and iteratively decides if a haplotype is matched based
on a translation table (in a gene)
"""

from . import logging


class MatchCounter:
    """
    Tests for each star allele
    """

    def __init__(self, haplotype, genotypes, chromosome, config):
        self.hap_positions = haplotype
        self.config = config
        self.chromosome = chromosome
        self.genotypes = genotypes
        self.missingsum = 0
        self.missings = []
        self.missingness = 0
        self.homsum = 0
        self.hetsum = 0
        self.matches = 0
        self.positions_matching = {}
        self.positions_matching_copy = {}
        self.varIDs_matched = []
        self.het_positions = []
        self.hom_positions = []
        self.position_counter = 0
        self.positions_not_matching = []
        self._iterate()
        if self.missingsum > 0:
            logging.warning(f"Missing markers for {self.hap_positions.name}: {', '.join(self.missings)}")

    def __str__(self):
        return self.hap_positions.name

    def __repr__(self):
        return self.hap_positions.name

    def update_match(self, allele):
        self.positions_matching_copy = self.positions_matching.copy()
        for key, item in allele.positions_matching.items():
            try:
                self.positions_matching_copy[key] -= 1
            except KeyError:
                pass

    def _iterate(self):
        for position in self.hap_positions.vars:
            try:
                subject = self.genotypes[f"{self.chromosome}:{position['Variant Start']}"]
            except KeyError:
                self.missings.append(position["Variant Start"])
                self.missingsum += 1
                continue
            alleles = subject["alleles"]
            if position["Type"] == "insertion":
                new_alleles = []
                for allele in subject["alleles"]:
                    if len(allele) == 1:
                        new_alleles.append("-")
                    else:
                        new_alleles.append(allele[-1])
                alleles = new_alleles
            self._test_alleles(alleles, position["Variant Allele"], position["Variant Start"])
    
    def _test_alleles(self, user_alleles, temp_alt, temp_pos):
        if temp_alt in self.config.IUPAC_CODES.keys():
            temp_alt = self.config.IUPAC_CODES[temp_alt]
        else:
            # put this in tuple so the 'in' test doesn't match part of string for in/del
            temp_alt = (temp_alt,)
        matches = sum([user_alleles.count(alt) for alt in temp_alt])
        if matches > 0:
            self._position_used_add(temp_pos, temp_alt, matches)
        else:
            self._position_not_used_add(temp_pos, temp_alt)

    def _position_used_add(self, position, nucleotide, matches):
        self.positions_matching[f"{position}-{nucleotide}"] = matches
        self.positions_matching_copy[f"{position}-{nucleotide}"] = matches

    def _position_not_used_add(self, position, nucleotide):
        self.positions_not_matching.append(f"{position}-{nucleotide}")

    def match(self):
        try:
            return True if len(self.positions_not_matching) == 0 and min([matches for position, matches in self.positions_matching_copy.items()]) > 0 else False
        except ValueError:
            return False
