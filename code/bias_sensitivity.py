# import externals libraries
import os
import numpy as np
import pandas as pd
import xarray as xr

import logging
log = logging.getLogger('bias-sensitivity')

# import the needed OGGM modules
from oggm import cfg, utils, workflow
from oggm.core import gis, climate, flowline
import oggm_vas as vascaling

import sys
sys.path.insert(0, os.getcwd())
from equilibrium_run import equilibrium_run_vas

# start logger with OGGM settings
cfg.set_logging_config()


equilibrium_run_vas(rgi_ids=['RGI60-11.00897'], use_random_mb=False,
                    path='../data/no-bias-sensitivity.nc')

