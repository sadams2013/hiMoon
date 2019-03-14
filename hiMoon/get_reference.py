from itertools import groupby


class GetReference(object):
    """
    Tools for getting reference positions from a fasta file
    """

    def __init__(self, fasta):
        self.fasta = fasta

    def chrom_contigs(self, chrom):
        fasta = (x[1] for x in groupby(
            open(self.fasta), lambda line: line[0] == ">"))
        fasta_data = {}
        i = 1
        for header in fasta:
            header = header.__next__()[1:].strip()
            seq = "".join(s.strip() for s in fasta.__next__())
            chrnum = header.split(":")[0]
            min = header.split(":")[1].split("-")[0]
            max = header.split(":")[1].split("-")[1]
            if str(chrnum[3:]) == str(chrom):
                fasta_data[i] = {"min": min, "max": max, "seq": seq}
            i += 1
        return fasta_data
