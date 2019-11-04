# Run the mass balance calibration

# Python imports
import json
import os

# Libs
import numpy as np

# Locals
import oggm
from oggm import cfg, utils, tasks, workflow
from oggm.workflow import execute_entity_task
from oggm.core import climate
from oggm.core.massbalance import (ConstantMassBalance, PastMassBalance,
                                   MultipleFlowlineMassBalance)


# RGI Version
rgi_version = '61'

# Since I'm in the Alps, I'll use HistAlp as baseline climate
baseline = 'HISTALP'

# Initialize OGGM and set up the run parameters
cfg.initialize()

# Local paths (where to write the OGGM run output)
dirname = 'OGGM_ref_mb_{}_RGIV{}_OGGM{}'.format(baseline, rgi_version,
                                                oggm.__version__)
WORKING_DIR = utils.gettempdir(dirname, home=True)
utils.mkdir(WORKING_DIR, reset=True)
cfg.PATHS['working_dir'] = WORKING_DIR


# The following code block alters certain parameters from
# the default config file.

# We are running the calibration ourselves
cfg.PARAMS['run_mb_calibration'] = True

# We are using which baseline data?
cfg.PARAMS['baseline_climate'] = baseline

# No need for intersects since this has an effect on the inversion only
cfg.PARAMS['use_intersects'] = False

# Use multiprocessing?
cfg.PARAMS['use_multiprocessing'] = False

# Set to True for operational runs
cfg.PARAMS['continue_on_error'] = False

if baseline == 'HISTALP':
    # Other params: see https://oggm.org/2018/08/10/histalp-parameters/
    # TODO: do I have to calibrate those parameters again, since I'm using
    # a different mass balance model?!
    cfg.PARAMS['baseline_y0'] = 1850
    cfg.PARAMS['prcp_scaling_factor'] = 1.75
    cfg.PARAMS['temp_melt'] = -1.75


# The next step is to get all the reference glaciers,
# i.e. glaciers with mass balance measurements.


# Get the reference glacier ids (they are different for each RGI version)
rgi_dir = utils.get_rgi_dir(version=rgi_version)
df, _ = utils.get_wgms_files()
rids = df['RGI{}0_ID'.format(rgi_version[0])]


# We can't do Antarctica
rids = [rid for rid in rids if not ('-19.' in rid)]

# For HISTALP only RGI reg 11
if baseline == 'HISTALP':
    rids = [rid for rid in rids if '-11.' in rid]

# Make a new dataframe with those (this takes a while)
print('Reading the RGI shapefiles...')
rgidf = utils.get_rgi_glacier_entities(rids, version=rgi_version)
print('For RGIV{} we have {} candidate reference '
      'glaciers.'.format(rgi_version, len(rgidf)))

# save those reference glaciers in a separate DataFrame
rgidf_alps = rgidf.copy()

# initialize the glacier regions
gdirs = workflow.init_glacier_regions(rgidf_alps)