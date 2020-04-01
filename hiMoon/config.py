import configparser
from . import logging

class ConfigData:

    def __init__(self, config_path):
        self.update_config = False
        # Initialize configparser
        self.config = configparser.ConfigParser()

        # Read config file
        self.config.read(config_path)

        try:
            self.CHROMOSOME_ACCESSIONS = {
                accession.upper(): chromosome for accession, chromosome in self.config["CHROMOSOME ACCESSIONS"].items()
                }
        except KeyError:
            logging.warning("Chromosome accessions not found in config, defaulting to GRCh38.")
            self.update_config = True
            # Defaults to GRCh38
            self.CHROMOSOME_ACCESSIONS = {
                "NC_000001.11": "1",
                "NC_000002.12": "2",
                "NC_000003.12": "3",
                "NC_000004.12": "4",
                "NC_000005.10": "5",
                "NC_000006.12": "6",
                "NC_000007.14": "7",
                "NC_000008.11": "8",
                "NC_000009.12": "9",
                "NC_000010.11": "10",
                "NC_000011.10": "11",
                "NC_000012.12": "12",
                "NC_000013.11": "13",
                "NC_000014.9": "14",
                "NC_000015.10": "15",
                "NC_000016.10": "16",
                "NC_000017.11": "17",
                "NC_000018.10": "18",
                "NC_000019.10": "19",
                "NC_000020.11": "20",
                "NC_000021.9": "21",
                "NC_000022.11": "22",
                "NC_000023.11": "23",
                "NC_000024.10": "24"
            }
            self.config["CHROMOSOME ACCESSIONS"] = self.CHROMOSOME_ACCESSIONS

        # Get IUPAC codes from config, if not found - use default and update
        try:
            self.IUPAC_CODES = {
                code.upper(): list(nucleotides) for code, nucleotides in self.config["IUPAC CODES"].items()
                }
        except KeyError:
            self.update_config = True
            logging.warning(
                "IUPAC Codes not found in config file. Adding default values.")
            self.IUPAC_CODES = {
                'R': ['A', 'G'], 
                'Y': ['C', 'T'], 
                'S': ['G', 'C'], 
                'W': ['A', 'T'], 
                'K': ['G', 'T'], 
                'M': ['A', 'C'], 
                'N': ['A', 'C', 'T', 'G']
                }
            self.config["IUPAC CODES"] = {
                code: "".join(nucleotides) for code, nucleotides in self.IUPAC_CODES.items()
                }

        # Get variant query parameters, if not found - add defaults
        try:
            self.VARIANT_QUERY_PARAMETERS = self.config["VARIANT QUERY PARAMETERS"]
        except KeyError:
            self.update_config = True
            logging.warning(
                "Query parameters not found in config file. Adding default values.")
            self.VARIANT_QUERY_PARAMETERS = {
                "5p_offset": 1000,
                "3p_offset": 1000
            }
            self.config["VARIANT QUERY PARAMETERS"] = self.VARIANT_QUERY_PARAMETERS

    