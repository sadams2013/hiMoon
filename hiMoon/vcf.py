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
import csv

import numpy as np

from pysam import VariantFile

from .template import PATH

from . import LOGGING


class VarFile:
    def __init__(self, vcf_file: str, sample: str = None, vcf_file_index: str = None, config = None) -> None:
        """
        VarFile object, basically a wrapper for pysam VariantFile
        
        Args:
            vcf_file (str): path to VCF/VCF.GZ/BCF file (needs to be indexed)
        """
        self.vcf_file = VariantFile(vcf_file, index_filename = vcf_file_index)
        if sample:
            self.samples = [sample]
        else:
            self.samples = list(self.vcf_file.header.samples)
        self.vcf_file.subset_samples(self.samples)
    
    def _get_alleles(self, sample, var_type):
        """
        Get alleles if coded as GT, catch if coded as CN (for CNV)
        Will always assume that at least one allele is one copy (unless CN == 0)
        Not necessarily accurate, but plain count data is inherently unphased so the distiction
        between (for example) 1/3 vs 2/2 is arbitrary. 
        """
        alleles = sample.alleles
        if len(alleles) == 0 and var_type == "CNV":
            try:
                cn = sample["CN"]
                if cn == 0:
                    alleles = ("0", "0")
                else:
                    alleles = ("1", str(cn - 1)) #TODO this in untested
            except KeyError:
                alleles = None # This is probably going to cause an issue eventually
        return alleles


    def get_range(self, chrom: str, minloc: int, maxloc: int) -> dict:
        """
        Returns a range of variants for all samples in a VCF file
        
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
            var_type = "SID"
            try:
                var_type = position.info["SVTYPE"]
            except KeyError:
                pass
            positions_out[f"c{chrom}_{position.pos}_{var_type}"] = {
                sample: {
                    "alleles": self._get_alleles(position.samples[sample], var_type), "phased": position.samples[sample].phased, "ref": position.ref} for sample in self.samples}
        return positions_out
        
def get_alleles(gene: object, subjects: list) -> list:
    """
    Prep for the ref/alt columns in a VCF

    Args:
        gene (Gene): gene object
        subjects (list): list of subjects

    Returns:
        list: list of all possible alt alleles
    """
    ref = gene.reference
    alts = []
    for s in subjects:
        subject_haps = [h for h in s.called_haplotypes[str(gene)]["HAPS"][0]]
        alts += np.array(subject_haps).flatten().tolist()
    try:
        alts = list(filter((ref).__ne__, alts))
    except ValueError:
        pass
    if len(alts) == 0:
        alts = ["NON_REF"]
    return([ref] + list(set(alts)))

def get_dosage(haps: list, alleles: list) -> list:
    """
    Get index for GT fiels for a given sample

    Args:
        haps (list): possible haplotypes
        alleles (list): alleles

    Returns:
        list: sample dosages for each allele
    """
    return [alleles.index(s) for s in haps]

def get_samples(gene_name: str, subjects: list, alleles: list) -> list:
    """
    Generate the sample/format field for each sample in the multi-sample output VCF

    Args:
        gene_name (str): gene name
        subjects (list): list of subjects to be written
        alleles (list): subject alleles

    Returns:
        list: list of format field for each each sample
    """
    formats = []
    for s in subjects:
        calls = s.called_haplotypes[gene_name]["HAPS"]
        formats.append(
            {
                "GT": get_dosage(calls[0][0], alleles),
                "VA": calls[1][0],
                "HC": 1 / len(calls[0])
            }
        )
    return formats


def write_variant_file(directory: str, subjects: [], prefix: str, genes: list) -> None:
    """
    Write the output VCF file

    Args:
        directory (str): output directory
        subjects ([type]): list of subjects
        prefix (str): prefix for filename
        genes (list): list of gene objects
    """
    contigs = list(set([f"chr{gene.chromosome.strip('chr')}" for gene in genes]))
    template = VariantFile(PATH + "/template.vcf", "r")
    outfile = VariantFile(directory+f"/{prefix}.haplotypes.vcf", "w", header = template.header)
    for contig in contigs:
        outfile.header.add_line(f"##contig=<ID={contig},length=0>")
    for sub in subjects:
        outfile.header.add_sample(str(sub))
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
            info = {"VARTYPE": "HAP"},
            samples = get_samples(str(gene), subjects, alleles)
        )
        outfile.write(nr)
    outfile.close()

def write_flat_file(directory: str, subjects: [], prefix: str) -> None:
    """
    Writes the output flat file

    Args:
        directory (str): output directory
        subjects ([type]): list of subjects
        prefix (str): prefix for filename
    """
    lines = []
    for subject in subjects:
        for gene, haps in subject.called_haplotypes.items():
            for i in range(len(haps["HAPS"][0])):
                lines.append({
                    "SUBJECT": str(subject),
                    "GENE": gene,
                    "GENOTYPE": "/".join(haps["HAPS"][0][i]),
                    "VARIANTS": "|".join(haps["HAPS"][1][i]),
                    "CONFIDENCE": 1 / len(haps["HAPS"][0])
                })
    with open(directory + f"/{prefix}.haplotypes.tsv", "w") as flat_out:
        flat_file = csv.DictWriter(flat_out, ["SUBJECT", "GENE", "GENOTYPE", "VARIANTS", "CONFIDENCE"], delimiter = "\t")
        flat_file.writeheader()
        flat_file.writerows(lines)
