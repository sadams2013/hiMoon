"""
Defines the Subject class
This serves to coordinate multiple gene lookups
for the same subject.
"""

import sys

from .vcf import VCFParse
from .match_counter import MatchCounter
from .alignment import AlignmentData
from .gene import Gene
from .pick_alleles import AllelePicker
from .cnv import find_cnv, find_cnv_simple

from . import logging, BAM_GENES, CNV_REGIONS, CNV_PENALTY


class Subject:
    """
    Describes a subject object
    Subject have called haplotypes and haplotype match objects
    These contain verbose information about how a haplotype was called
    """

    def __init__(
        self, 
        prefix: str,
        genes: [], 
        vcf_file: str = None, 
        alignment_file: str = None, 
        copy_number_control: str = None):
        self.prefix = prefix
        self.called_haplotypes = {}
        self.cnv_regions = []
        self.genes = {}
        self.genotypes = {}
        self.bam = None
        self.vcf = None
        self.copy_number_control = None
        self.normalization_factors = {}
        if vcf_file:
            self.vcf = VCFParse(vcf_file = vcf_file)
        if alignment_file:
            self.bam = AlignmentData(alignment_file = alignment_file)
        if copy_number_control:
            self.copy_number_control = AlignmentData(alignment_file = copy_number_control)
        for gene in genes:
            genotypes = self.get_gene(gene)
            haplotypes = {}
            for name, haplotype in gene.haplotypes.items():
                haplotypes[name] = MatchCounter(haplotype, genotypes, gene.chromosome())
            self.genes[str(gene)] = haplotypes
            self.called_haplotypes[str(gene)] = AllelePicker(gene, haplotypes).pick_stars()
        if alignment_file and copy_number_control:
            self.cnv()
        elif alignment_file:
            self.cnv_simple()
            
    
    def cnv(self):
        cnv_regions = []
        for region, data in CNV_REGIONS.items():
                cnv_data = find_cnv(
                    case = self.bam.get_read_array(data[0], int(data[1]), int(data[2])),
                    control = self.copy_number_control.get_read_array(data[0], int(data[1]), int(data[2])),
                    penalty = CNV_PENALTY)
                ss = {}
                used = []
                for span in cnv_data:
                    if span[2] <= 0.6 or span[2] >= 1.4: 
                        start = span[0]
                        stop = span[1]
                        rat = span[2]
                        if start not in used and int(stop) > int(start):
                            ss[start] = (stop, rat)
                for start, span_data in ss.items():
                    cnv_regions.append({
                        "PREFIX": self.prefix,
                        "REGION": region,
                        "CHROMOSOME": data[0],
                        "SPAN_START": int(start),
                        "SPAN_STOP": int(span_data[0]),
                        "LENGTH": int(span_data[0]) - int(start),
                        "SPAN_RATIO": span_data[1]
                    })
        self.cnv_regions = self._simplify_cnv(cnv_regions)


    def cnv_simple(self):
        for region, data in CNV_REGIONS.items():
            try:
                span_ratio = find_cnv_simple(
                    case = self.bam.get_read_array(data[0], int(data[1]), int(data[2])),
                    control = self.bam.get_read_array(data[0], int(data[3]), int(data[4]))
                )
                self.cnv_regions.append({
                        "PREFIX": self.prefix,
                        "REGION": region + "_simple",
                        "CHROMOSOME": data[0],
                        "SPAN_START": data[1],
                        "SPAN_STOP": data[2],
                        "LENGTH": int(data[2]) - int(data[1]),
                        "SPAN_RATIO": span_ratio
                })
            except IndexError:
                continue

    
    def get_gene(self, gene: Gene) -> dict:
        """
        Wrapper for vcf parser and bam parser
        """
        if str(gene) in BAM_GENES and self.bam or not self.vcf:
            return self.bam.get_range(
                                    gene.chromosome(),
                                    gene.position_min(), 
                                    gene.position_max()
                                    )
        elif self.vcf:
            return self.vcf.get_range(
                                    gene.chromosome(),
                                    gene.position_min(), 
                                    gene.position_max()
                                    )

    def _simplify_cnv(self, cnvs: []) -> list:
        """
        Limit returned CNV regions to the largest with no overlap. 
        """
        new_cnvs = []
        for cnv in cnvs:
            keep = True
            for cnv_ in cnvs:
                if cnv_["CHROMOSOME"] == cnv["CHROMOSOME"]: 
                    if cnv_["SPAN_START"] <= cnv["SPAN_START"]:
                        if cnv_["SPAN_STOP"] > cnv["SPAN_STOP"]:
                            keep = False
                    if cnv_["SPAN_STOP"] >= cnv["SPAN_STOP"]:
                        if cnv_["SPAN_START"] < cnv["SPAN_START"]:
                            keep = False
            if keep:
                new_cnvs.append(cnv)
        return new_cnvs
