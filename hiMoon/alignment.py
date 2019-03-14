import sys
import math
import numpy as np

from pysam import AlignmentFile
from functools import reduce
from . import MAX_PILEUP_DEPTH, CALL_THRESHOLD, MIN_READS, MIN_Q, MIN_MAPQ

from .gene import Gene

class AlignmentData:

    def __init__(self, alignment_file: str, type = "b"):
        self.alignment_file = alignment_file
        self.alignment_object = AlignmentFile(alignment_file, f"r{type}")
        try:
            self.alignment_object.check_index()
        except AttributeError:
            print("SAM files must be converted to BAM and indexed!")
            sys.exit(1)
        except ValueError:
            print(
                "No index file found! Please index with 'samtools index <file.bam/cram>'")
            sys.exit(1)

    def get_range(self, chrom: str, start: int, stop: int) -> list:
        """Returns range based on chromosome and min/max locations"""
        positions_data = self.call_alleles(chrom, start, stop)
        return positions_data
    
    def call_alleles(
        self, 
        chrom: str, 
        start: int, 
        stop: int, 
        min_cov: int = MIN_READS, 
        cutoff_freq: float = CALL_THRESHOLD
        ) -> dict:
        positions_out = {}
        gene_pileups = self.get_pileup( chrom, start, stop)
        for pileupcolumn in gene_pileups:
            position = pileupcolumn.reference_pos
            positions_out[f"{chrom}:{position + 1}"] = self._call_from_reads(
                    pileupcolumn.get_query_sequences(add_indels=True), 
                    min_cov, 
                    cutoff_freq
                    )
        return positions_out


    def get_read_array(self, chromosome: str, start: int, stop: int) -> np.ndarray:
        starts = []
        stops = []
        for read in self._fetch_wrapper(chromosome, start, stop):
            read_list = str(read).split("\t")
            start = int(read_list[3])
            length_of_segment = int(read_list[8])
            starts.append(start)
            stops.append(start + length_of_segment)
        return np.array(starts, dtype=np.int32)

    def get_pileup(
        self, 
        chrom: str, 
        start: int, 
        stop: int,
        min_mapping_quality: int = MIN_MAPQ,
        min_base_quality: int = MIN_Q,
        max_depth: int = MAX_PILEUP_DEPTH) -> list:
        try:
            return self.alignment_object.pileup(
                contig = str(chrom), 
                start = start, 
                stop = stop, 
                truncate = True,
                max_depth = max_depth,
                min_mapping_quality = min_mapping_quality,
                min_base_quality = min_base_quality
        )
        except ValueError:
            return self.alignment_object.pileup(
                contig = f"chr{chrom}", 
                start = start, 
                stop = stop,
                truncate = True,
                max_depth = max_depth,
                min_mapping_quality = min_mapping_quality,
                min_base_quality = min_base_quality
                )

    def _fetch_wrapper(self, chrom, start, stop):
        try:
            return self.alignment_object.fetch(
                contig = str(chrom), 
                start = start, 
                stop = stop)
        except ValueError:
            return self.alignment_object.fetch(
                contig = f"chr{chrom}", 
                start = start, 
                stop = stop)


    def _call_from_reads(self, pileup_query: list, min_cov: float, cutoff_freq: float) -> dict:
        """
        Make allele calls based on pileup. Calculates the fraction of all reads that map to 
        each allele and returns estimated allele calls.
        """
        # Get fraction of all possible alleles. Converts all to upper case for simplicity.
        alleles = []
        for allele in set(pileup_query):
            if pileup_query.count(allele)/len(pileup_query) > cutoff_freq:
                if allele.upper() not in ["A", "C", "T", "G", "*"]:
                    if "+" in allele:
                        allele = "".join([nuc for nuc in allele.replace("+", "") if not nuc.isdigit()])
                    if "-" in allele:
                        allele = allele[0]
                allele = allele.replace("*", "-")
                alleles.append(allele.upper())
        # Add the other allele if homozygous
        if len(alleles) == 1:
            alleles += alleles[0]
        elif len(alleles) > 2:
            alleles = list(set(alleles))
        return {
                "alleles": alleles,
                "phased": False,
                "ref": "",
                "ID": "NA",
                "source": "alignment",
                "relative_coverage": None
                }
