#    Copyright 2021 Solomon M. Adams

#    Licensed under the Apache License, Version 2.0 (the "License");
#    you may not use this file except in compliance with the License.
#    You may obtain a copy of the License at

#        http://www.apache.org/licenses/LICENSE-2.0

#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#    See the License for the specific language governing permissions and
#    limitations under the License.

__author__ = "Solomon M. Adams, PharmD, PhD"
__copyright__ = "Copyright 2021, Solomon M. Adams"
__license__ = "Apache v2.0"
__version__ = "0.13.2"
__maintainer__ = "Solomon M. Adams, PharmD, PhD"
__email__ = "sadams2013@gmail.com"

import logging as LOGGING

from .config import ConfigData

SPECIAL_CHROM = {
    "23": "X",
    "24": "Y"
}

# Initialize logger
LOGGING.getLogger().setLevel(LOGGING.WARNING)

def set_logging_info():
    LOGGING.getLogger().setLevel(LOGGING.INFO)

def get_config(config_path: str = None):
    """
    Set a custom config based on a path
    Args:
        config_path ([type]): path/to/config/file.ini
    """
    if config_path:
        return(ConfigData(config_path))
    else:
        return(ConfigData())


