import configparser

from . import LOGGING

GRCH38_ACCESSIONS = {
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

IUPAC_CODES = {
                'R': ['A', 'G'], 
                'Y': ['C', 'T'], 
                'S': ['G', 'C'], 
                'W': ['A', 'T'], 
                'K': ['G', 'T'], 
                'M': ['A', 'C'], 
                'N': ['A', 'C', 'T', 'G'],
                'A': ['A'],
                'C': ['C'],
                'T': ['T'],
                'G': ['G']                
            }

class ConfigData:

    def __init__(self, config_path: str = "config.ini", write_config: bool = False) -> None:
        """
        Create a new config class that is used throughout to set various parameters

        Args:
            config_path (str): path to a configuration file - can be a "dummy" value. 
            write_config (bool, optional): Write a new configuration file to config_path?. Defaults to False.
        """
        self.config = configparser.ConfigParser()
        self.config.read(config_path)
        self.chromosome_accessions()
        self.iupac_codes()
        self.variant_query_params()
        if write_config:
            self.write_config(config_path)
    
    def chromosome_accessions(self) -> None:
        """
        Set chromosome accessions
        """
        try:
            self.CHROMOSOME_ACCESSIONS = {accession.upper(): chromosome for accession, chromosome in self.config["CHROMOSOME ACCESSIONS"].items()}
            print(self.CHROMOSOME_ACCESSIONS)
        except KeyError:
            LOGGING.info("Chromosome accessions defaulting to GRCh38")
            self.CHROMOSOME_ACCESSIONS = GRCH38_ACCESSIONS
            self.config["CHROMOSOME ACCESSIONS"] = self.CHROMOSOME_ACCESSIONS
    
    def iupac_codes(self) -> None:
        """
        Set IUPAC nucleotide codes
        """
        # Get IUPAC codes from config, if not found - use default and update
        try:
            self.IUPAC_CODES = {code.upper(): list(nucleotides) for code, nucleotides in self.config["IUPAC CODES"].items()}
        except KeyError:
            self.update_config = True
            LOGGING.info("IUPAC Codes set to default values.")
            self.IUPAC_CODES = IUPAC_CODES
            self.config["IUPAC CODES"] = {code: "".join(nucleotides) for code, nucleotides in self.IUPAC_CODES.items()}
    
    def variant_query_params(self) -> None:
        """
        Set variant query parameters
        """
        # Get variant query parameters, if not found - add defaults
        try:
            self.VARIANT_QUERY_PARAMETERS = self.config["VARIANT QUERY PARAMETERS"]
        except KeyError:
            self.update_config = True
            LOGGING.info("Using a 1kb 5' and 3' offset.")
            self.VARIANT_QUERY_PARAMETERS = {
                "5p_offset": 1000,
                "3p_offset": 1000
            }
            self.config["VARIANT QUERY PARAMETERS"] = self.VARIANT_QUERY_PARAMETERS
    
    def write_config(self, config_path: str) -> None:
        """
        Write default/modified parameters to file at config_path

        Args:
            config_path (str): path to configuration file
        """
        with open(config_path, "w") as configfile:	
            self.config.write(configfile)
    