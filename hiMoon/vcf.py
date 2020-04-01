"""
Reads a VCF (gz) file, then allows the user to lookup by range.
File must be bgzipped and tabixed
"""

from pysam import VariantFile

from . import logging


class VarFile:
    def __init__(self, vcf_file: str):
        """
        Create a new VCF file object.
        """
        self.vcf_file = VariantFile(vcf_file)
        self.samples = list(self.vcf_file.header.samples)
    
    def get_range(self, chrom: str, minloc: int, maxloc: int) -> dict:
        """Returns range based on chromosome and min/max locations"""
        positions_out = {}
        try:
            positions = self.vcf_file.fetch(str(chrom), minloc, maxloc)
        except ValueError:
            positions = self.vcf_file.fetch(f"chr{chrom}", minloc, maxloc)
        for position in positions:
            chrom = position.chrom.strip("chr")
            positions_out[f"c{chrom}_{position.pos}"] = {
                sample: {
                    "alleles": tuple(position.samples[sample].alleles), "phased": position.samples[sample].phased, "ref": position.ref} for sample in self.samples}
        return positions_out
