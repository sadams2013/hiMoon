"""
Reads a VCF (gz) file, then allows the user to lookup by range.
File must be bgzipped and tabixed
"""

from pysam import VariantFile

from . import logging


class VarFile:
    def __init__(self, vcf_file: str, sample: str = None) -> None:
        """VarFile object, basically a wrapper for pysam VariantFile
        
        Args:
            vcf_file (str): path to VCF/VCF.GZ/BCF file (needs to be indexed)
        """
        self.vcf_file = VariantFile(vcf_file)
        if sample:
            self.samples = [sample]
        else:
            self.samples = list(self.vcf_file.header.samples)
        self.vcf_file.subset_samples(self.samples)
    
    def get_range(self, chrom: str, minloc: int, maxloc: int) -> dict:
        """Returns a range of variants for all samples in a VCF file
        
        Args:
            chrom (str): chromosome
            minloc (int): starting position
            maxloc (int): ending position
        
        Returns:
            dict: variants with a common ID schema that is matched by other methods
        """
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
