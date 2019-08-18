"""
Reads a VCF (gz) file, then allows the user to lookup by range.
File must be bgzipped and tabixed
"""

from pysam import VariantFile

from . import logging


class VCFParse:
    def __init__(self, vcf_file: str, index = None):
        """
        Create a new VCF file object.
        """
        self.vcf_file = VariantFile(vcf_file, index_filename = index)

    def get_range(self, chrom: str, minloc: int, maxloc: int) -> dict:
        """Returns range based on chromosome and min/max locations"""
        positions_out = {}
        try:
            positions = self.vcf_file.fetch(
                str(chrom), minloc, maxloc, reopen=True)
        except ValueError:
            positions = self.vcf_file.fetch(
                f"chr{chrom}", minloc, maxloc, reopen=True)
        for position in positions:
            sample = position.samples[0]
            chrom = position.chrom.strip("chr")
            positions_out[f"{chrom}:{position.pos}"] = {"alleles": tuple(sample.alleles),
                                                        "phased": sample.phased, 
                                                        "ID": position.id, 
                                                        "ref": position.ref,
                                                        "source": "vcf",
                                                        "relative_coverage": None,
                                                        "structural": False}
        return positions_out
