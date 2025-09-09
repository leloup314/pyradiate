#!/usr/bin/env python

import os
import logging
from pathlib import Path

#from pyradiate.tools import check_install


# Version attribute that is picked dynaimcally by setuptools
__version__ =  "0.0.1"


# Package paths
package_path = Path(__file__).parent
data_path = Path.home() / Path(".local/share/pyradiate")
config_path = Path.home() / Path(".config/pyradiate")    
ensdf_path = Path.home() / Path(".local/share/pyradiate/ensdf")

# Package files
ensdf_file = data_path / "ensdf.zip"
ensdf_config = config_path / "ensdf.yaml"

# Configure path for nudel
os.environ["ENSDF_PATH"] = str(ensdf_path)

# Check if required paths exists; if not, create
if not data_path.exists():
    data_path.mkdir()

if not config_path.exists():
    config_path.mkdir()

if not ensdf_path.exists():
    ensdf_path.mkdir()


# Create logger
logging.basicConfig(level=logging.DEBUG,
                    format="%(asctime)s [%(name)-17s] - %(levelname)-7s %(message)s")
logger = logging.getLogger("pyRADIATE")