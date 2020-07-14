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

import configparser

from . import LOGGING

BASE_ACCESSION = {
    ".": "NA"
}

GRCH37_ACCESSIONS = {
                "NC_000001.10": "1",
                "NC_000002.11": "2",
                "NC_000003.11": "3",
                "NC_000004.11": "4",
                "NC_000005.9": "5",
                "NC_000006.11": "6",
                "NC_000007.13": "7",
                "NC_000008.10": "8",
                "NC_000009.11": "9",
                "NC_000010.10": "10",
                "NC_000011.9": "11",
                "NC_000012.11": "12",
                "NC_000013.10": "13",
                "NC_000014.8": "14",
                "NC_000015.9": "15",
                "NC_000016.9": "16",
                "NC_000017.10": "17",
                "NC_000018.9": "18",
                "NC_000019.9": "19",
                "NC_000020.10": "20",
                "NC_000021.8": "21",
                "NC_000022.10": "22",
                "NC_000023.10": "23",
                "NC_000024.9": "24"
            }

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

    def __init__(self, config_path: str = "config.ini") -> None:
        """
        Create a new config class that is used throughout to set various parameters

        Args:
            config_path (str): path to a configuration file - can be a "dummy" value. 
            write_config (bool, optional): Write a new configuration file to config_path?. Defaults to False.
        """
        self.config = configparser.ConfigParser()
        self.config.read(config_path)
        self._chromosome_accessions()
        self._iupac_codes()
        self._variant_query_params()
        self._missing_params()
        self._lp_params()
    
    def _chromosome_accessions(self) -> None:
        """
        Set chromosome accessions
        """
        try:
            self.CHROMOSOME_ACCESSIONS = {accession.upper(): chromosome for accession, chromosome in self.config["CHROMOSOME ACCESSIONS"].items()}
        except KeyError:
            LOGGING.info("Chromosome accessions defaulting to GRCh38")
            self.CHROMOSOME_ACCESSIONS = {**GRCH38_ACCESSIONS, **GRCH37_ACCESSIONS}
            self.config["CHROMOSOME ACCESSIONS"] = self.CHROMOSOME_ACCESSIONS
        self.CHROMOSOME_ACCESSIONS = {**self.CHROMOSOME_ACCESSIONS, **BASE_ACCESSION}
    
    def _iupac_codes(self) -> None:
        """
        Set IUPAC nucleotide codes
        """
        # Get IUPAC codes from config, if not found - use default and update
        try:
            self.IUPAC_CODES = {code.upper(): list(nucleotides) for code, nucleotides in self.config["IUPAC CODES"].items()}
        except KeyError:
            LOGGING.info("IUPAC Codes set to default values.")
            self.IUPAC_CODES = IUPAC_CODES
            self.config["IUPAC CODES"] = {code: "".join(nucleotides) for code, nucleotides in self.IUPAC_CODES.items()}
    
    def _variant_query_params(self) -> None:
        """
        Set variant query parameters
        """
        # Get variant query parameters, if not found - add defaults
        try:
            self.VARIANT_QUERY_PARAMETERS = self.config["VARIANT QUERY PARAMETERS"]
        except KeyError:
            LOGGING.info("Using a 1kb 5' and 3' offset.")
            self.VARIANT_QUERY_PARAMETERS = {
                "5p_offset": 1000,
                "3p_offset": 1000
            }
            self.config["VARIANT QUERY PARAMETERS"] = self.VARIANT_QUERY_PARAMETERS
    
    def _missing_params(self) -> None:
        try:
            self.MISSING_DATA_PARAMETERS = self.config["MISSING DATA PARAMETERS"]
        except KeyError:
            LOGGING.info("Missing variants will be assumed to match reference")
            self.MISSING_DATA_PARAMETERS = {
                "missing_variants": 99
            }
            self.config["MISSING DATA PARAMETERS"] = self.MISSING_DATA_PARAMETERS
    
    def _lp_params(self) -> None:
        try:
            self.LP_PARAMS = self.config["LINEAR PROGRAM PARAMETERS"]
        except KeyError:
            self.LP_PARAMS = {
                "optimal_decay": 0,
                "max_haps": 2
            }
            self.config["LINEAR PROGRAM PARAMETERS"] = self.LP_PARAMS
    
    def write_config(self, config_path: str) -> None:
        """
        Write default/modified parameters to file at config_path

        Args:
            config_path (str): path to configuration file
        """
        with open(config_path, "w") as configfile:	
            self.config.write(configfile)
    