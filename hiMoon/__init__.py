__author__ = "Solomon M. Adams, PharmD, PhD"
__copyright__ = "Copyright 2020, Solomon M. Adams"
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Solomon M. Adams, PharmD, PhD"
__email__ = "sadams07@su.edu"

import logging as LOGGING

from .config import ConfigData

# Initialize logger
LOGGING.getLogger().setLevel(LOGGING.WARNING)

# Initialize config
CONFIG = ConfigData()

def set_logging_info():
    LOGGING.getLogger().setLevel(LOGGING.INFO)

def set_config(config_path):
    """
    Set a custom config based on a path
    Args:
        config_path ([type]): path/to/config/file.ini
    """
    CONFIG = ConfigData(config_path)
