import csv

"""
Take reference FASTA and impute reference positions that return like vcf_parser objects
"""


class FastaPositions(object):

    def __init__(self, fasta):
        self.contigs = []
        self.fasta_data = self._parse_fasta(fasta)

    def get_ref(self, chrom, position):
        try:
            return self.fasta_data[f"{chrom}:{position}"]
        except KeyError:
            return {"allele1": "", "allele2": "", "phased": False, "source": "fasta", "ref": "N"}

    def _parse_fasta(self, fasta_path):
        """
        :param fasta_path:path to fasta file
        :return: dictionary containing position, chromosome, and ref allele
        """
        outs = {}
        with open(fasta_path, "r") as fasta_file:
            fasta_lines = iter(csv.reader(fasta_file))
            for line in fasta_lines:
                chrom = line[0].split(":")[0].strip(">chr")
                position_start = int(line[0].split(":")[1].split("-")[0]) + 1
                position_end = int(line[0].split(":")[1].split("-")[1]) + 1
                self.contigs.append(
                    {"start": position_start, "end": position_end, "chrom": chrom})
                positions = list(range(position_start, position_end))
                alleles = next(fasta_lines)[0]
                for position in positions:
                    outs[f"{chrom}:{position}"] = {"ID": "", "allele1": alleles[position - position_start],
                                                   "allele2": alleles[position - position_start],
                                                   "phased": True, "ref": alleles[position - position_start], "source": "fasta"}
        return outs
