__author__ = "Solomon M. Adams, PharmD, PhD"
__copyright__ = "Copyright 2019, Solomon M. Adams"
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Solomon M. Adams, PharmD, PhD"
__email__ = "sadams07@su.edu"

import configparser
import logging

config_file = "config.ini"

# Run the config update, defaults to false
update_config = False

# Initialize logger
logging.getLogger().setLevel(logging.WARNING)

# Initialize configparser
config = configparser.ConfigParser()

# Read config file
config.read(config_file)

try:
    CHROMOSOME_ACCESSIONS = {
        accession.upper(): chromosome for accession, chromosome in config["CHROMOSOME ACCESSIONS"].items()
        }
except KeyError:
    logging.warning("Chromosome accessions not found in config, defaulting to GRCh38.")
    update_config = True
    # Defaults to GRCh38
    CHROMOSOME_ACCESSIONS = {
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
    config["CHROMOSOME ACCESSIONS"] = CHROMOSOME_ACCESSIONS

# Get IUPAC codes from config, if not found - use default and update
try:
    IUPAC_CODES = {
        code.upper(): list(nucleotides) for code, nucleotides in config["IUPAC CODES"].items()
        }
except KeyError:
    update_config = True
    logging.warning(
        "IUPAC Codes not found in config file. Adding default values.")
    IUPAC_CODES = {
        'R': ['A', 'G'], 
        'Y': ['C', 'T'], 
        'S': ['G', 'C'], 
        'W': ['A', 'T'], 
        'K': ['G', 'T'], 
        'M': ['A', 'C'], 
        'N': ['A', 'C', 'T', 'G']
        }
    config["IUPAC CODES"] = {
        code: "".join(nucleotides) for code, nucleotides in IUPAC_CODES.items()
        }

# Get variant query parameters, if not found - add defaults
try:
    VARIANT_QUERY_PARAMETERS = config["VARIANT QUERY PARAMETERS"]
except KeyError:
    update_config = True
    logging.warning(
        "Query parameters not found in config file. Adding default values.")
    VARIANT_QUERY_PARAMETERS = {
        "5p_offset": 1000,
        "3p_offset": 1000
    }
    config["VARIANT QUERY PARAMETERS"] = VARIANT_QUERY_PARAMETERS

# Look for bam genes in config, otherwise set to empty
try: 
    BAM_GENES = config["BAM GENES"]["bam_file_genes"].split(",")
except KeyError:
    logging.info(
        "No genes included in bam genes."
    )
    BAM_GENES = None

if update_config:
    with open(config_file, "w") as configfile:
        config.write(configfile)


## All following factors DO NOT write to the config.ini that is set by this file. If you need to know how to format those
## See the default config file in the repository. 

# Get maximum pileup depth
try: 
    MAX_PILEUP_DEPTH = int(config["BAM GENES"]["max_pileup_depth"])
    CALL_THRESHOLD = float(config["BAM GENES"]["call_threshold"])
    MIN_READS = int(config["BAM GENES"]["min_reads"])
    MIN_Q = int(config["BAM GENES"]["min_q"])
    MIN_MAPQ = int(config["BAM GENES"]["min_mapq"])
except KeyError:
    logging.info(
        "Missing one or more variant calling parameters, setting all to default."
    )
    MAX_PILEUP_DEPTH = 8000
    CALL_THRESHOLD = 0.2
    MIN_READS = 30
    MIN_Q = 30
    MIN_MAPQ = 30

# Look for CNV information in config. Does not write a default. 

try:
    CNV_REGIONS = {name: tuple(data.split(",")) for name, data in config["CNV REGIONS"].items()}
except KeyError:
    CNV_REGIONS = None

try:
    CNV_PENALTY = int(config["CNV PARAMETERS"]["penalty"])
except KeyError:
    CNV_PENALTY = 5
