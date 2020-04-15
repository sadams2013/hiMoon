"""
Reads a VCF (gz) file, then allows the user to lookup by range.
File must be bgzipped and tabixed
"""

from pysam import VariantFile

from .template import PATH

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

def get_alleles(gene, subjects):
    ref = gene.reference
    alts = []
    for s in subjects:
        subject_haps = [h for h in s.called_haplotypes[str(gene)]["HAPS"][1]]
        alts += subject_haps
    try:
        alts.remove(ref)
    except ValueError:
        pass
    return([ref] + list(set(alts)))

def get_dosage(haps, alleles):
    return [alleles.index(s) for s in haps]


def get_samples(gene_name, subjects, alleles):
    formats = []
    for s in subjects:
        formats.append(
            {
                "GT": get_dosage(s.called_haplotypes[gene_name]["HAPS"][1], alleles),
                "VA": s.called_haplotypes[gene_name]["HAPS"][2]
            }
        )
    return formats


def write_variant_file(directory: str, subjects: [], prefix, genes):
    contigs = list(set([f"chr{gene.chromosome.strip('chr')}" for gene in genes]))
    template = VariantFile(PATH + "/template.vcf", "r")
    outfile = VariantFile(directory+f"/{prefix}.haplotypes.vcf", "w", header = template.header)
    for contig in contigs:
        outfile.header.add_line(f"##contig=<ID={contig},length=0>")
    for sub in subjects:
        outfile.header.add_sample(str(sub))
    records = {}
    for gene in genes:
        alleles = get_alleles(gene, subjects)
        nr =outfile.new_record(
            contig = f"chr{gene.chromosome}",
            start = gene.min,
            stop = gene.max,
            alleles = [f'<{a.replace(str(gene), "").replace("(star)", "*")}>' for a in alleles],
            id = f"{str(gene)}_pgx",
            qual = None,
            filter = None,
            info = None,
            samples = get_samples(str(gene), subjects, alleles)
        )
        outfile.write(nr)
    outfile.close()
