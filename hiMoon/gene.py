"""
Defines the class Gene and Haplotype
Gene contains methods useful for the gene
Haplotype contains methods for gene specific haplotypes
"""

import csv

from . import VARIANT_QUERY_PARAMETERS, BAM_GENES, CHROMOSOME_ACCESSIONS, config, logging


class Gene:
    """
    Describes a gene object and provides clean library data.
    This contains pertinent information about the gene
    and its haplotypes
    """

    def __init__(self, translation_table):
        self.haplotypes = {}
        self.gene = None
        self.accession = None
        self.ref = "Ref"
        with open(translation_table, 'rt') as trans_file:
            self.version = trans_file.readline().strip("#version=\n\t")
            haplotype_file_lines = csv.DictReader(trans_file, delimiter="\t")
            for hap_line in haplotype_file_lines:
                self.gene = hap_line["Gene"]
                try:
                    if self.haplotypes[
                        hap_line["Haplotype Name"]
                        ].add_var(hap_line):
                        self.ref = hap_line["Haplotype Name"]
                    else:
                        self.accession = hap_line["ReferenceSequence"]
                except KeyError:
                    self.haplotypes[
                        hap_line["Haplotype Name"]
                        ] = Haplotype()
                    if self.haplotypes[
                        hap_line["Haplotype Name"]
                        ].add_var(hap_line):
                        self.ref = hap_line["Haplotype Name"]
                    else:
                        self.accession = hap_line["ReferenceSequence"]

    def __str__(self):
        return self.gene

    def __repr__(self):
        return self.gene

    def chromosome(self) -> int:
        return int(CHROMOSOME_ACCESSIONS[
            self.accession
        ])

    def position_min(self, offset=VARIANT_QUERY_PARAMETERS["5p_offset"]) -> int:
        return min([hap.get_min() for name, hap in self.haplotypes.items()]) - int(offset)

    def position_max(self, offset=VARIANT_QUERY_PARAMETERS["3p_offset"]) -> int:
        return max([hap.get_max() for name, hap in self.haplotypes.items()]) + int(offset)


class Haplotype:
    """
    Haplotype objects - sorts and stores characteristics of markers
    This class is only instantiated from within the gene class
    """

    def __init__(self):
        self.name = None
        self.is_ref = False
        self.ref_row = None
        self.vars = []

    def add_var(self, var_data: dict) -> bool:
        """
        Takes a row from the 
        """
        self.name = var_data["Haplotype Name"]
        if var_data["ReferenceSequence"] == "REFERENCE":
            self.is_ref = True
            self.ref_row = var_data
        else:
            self.vars.append(var_data)
        return self.is_ref

    def get_min(self) -> int:
        """Returns the minimum position needed for this haplotype"""
        try:
            return min(
                [int(variant["Variant Start"]) for variant in self.vars]
                )
        except ValueError:
            return 100000000

    def get_max(self) -> int:
        """Returns the maximum position needed for this haplotype"""
        try:
            return max(
                [int(variant["Variant Stop"]) for variant in self.vars]
                )
        except ValueError:
            return 0
