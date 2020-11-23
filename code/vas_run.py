""" Standard VAS run using the example of Hintereisferner.
Can be used as template.
"""
# Built ins
import os
import json
import logging
import datetime
from time import gmtime, strftime

# External libs
import numpy as np
import pandas as pd
import xarray as xr
import netCDF4
from scipy.optimize import minimize_scalar

# import OGGM modules
import oggm_vas as vascaling
import oggm
from oggm import cfg, utils, workflow, tasks
from oggm.shop import rgitopo
from oggm.cfg import SEC_IN_YEAR, SEC_IN_MONTH

from oggm import __version__

from oggm.utils import floatyear_to_date, ncDataset
from oggm.exceptions import InvalidParamsError, MassBalanceCalibrationError

from oggm.core import climate, gis
from oggm.core.massbalance import MassBalanceModel

# Module logger
log = logging.getLogger(__name__)

# define the glaciers you want to use by their RGI ID
rgi_ids = ['RGI60-11.00897']

# OGGM preprocessing
# ------------------
# load default parameter file
vascaling.initialize()
log.info('Starting preprocessiog')

# compute RGI region and version from RGI IDs
# assuming they all are all the same
rgi_region = (rgi_ids[0].split('-')[-1]).split('.')[0]
rgi_version = (rgi_ids[0].split('-')[0])[-2:-1]

# specify working directory and output directory
working_dir = oggm.utils.get_temp_dir('vas_run')
log.info(f'Working directory for this run: {working_dir}')
# output_dir = os.path.abspath('./vas_run_output')
output_dir = os.path.abspath('../data/vas_run_output')
# create working directory
utils.mkdir(working_dir)
utils.mkdir(output_dir)
# set path to working directory
cfg.PATHS['working_dir'] = working_dir
# set RGI version and region
cfg.PARAMS['rgi_version'] = rgi_version
# define how many grid points to use around the glacier,
# if you expect the glacier to grow large use a larger border
cfg.PARAMS['border'] = 20
# we use HistAlp climate data
cfg.PARAMS['baseline_climate'] = 'HISTALP'
# set the mb hyper parameters accordingly
cfg.PARAMS['prcp_scaling_factor'] = 2.5
cfg.PARAMS['temp_melt'] = -0.5
cfg.PARAMS['run_mb_calibration'] = False

# the bias is defined to be zero during the calibration process,
# which is why we don't use it here to reproduce the results
cfg.PARAMS['use_bias_for_run'] = True

# read RGI entry for the glaciers as DataFrame
# containing the outline area as shapefile
rgidf = utils.get_rgi_glacier_entities(rgi_ids)

# get and set path to intersect shapefile
intersects_db = utils.get_rgi_intersects_region_file(region=rgi_region)
cfg.set_intersects_db(intersects_db)

# sort by area for more efficient parallel computing
rgidf = rgidf.sort_values('Area', ascending=False)
cfg.PARAMS['use_multiprocessing'] = True
# operational run, all glaciers should run
cfg.PARAMS['continue_on_error'] = False

# initialize the GlacierDirectory
gdirs = workflow.init_glacier_directories(rgi_ids)

# define the local grid and glacier mask
# TODO: solv problem with NASADEM as source
workflow.execute_entity_task(tasks.define_glacier_region, gdirs, source='SRTM')
workflow.execute_entity_task(gis.glacier_masks, gdirs)
# process the given climate file
workflow.execute_entity_task(climate.process_climate_data, gdirs)
# compute local t* and the corresponding mu*
workflow.execute_entity_task(vascaling.local_t_star, gdirs)

# --------------------
#  SCALING MODEL
# --------------------

# compute local t* and the corresponding mu*
log.info('Initialize the model')

gdir = gdirs[0]

# instance the mass balance models
mb_mod = vascaling.VAScalingMassBalance(gdir)

# get reference area
a0 = gdir.rgi_area_m2
# get reference year
y0 = int(rgidf.BgnDate.values[0][:4])
# get min and max glacier surface elevation
h0, h1 = vascaling.get_min_max_elevation(gdir)

# set up the glacier model with the values of 2003
hef_ref = vascaling.VAScalingModel(year_0=1851, area_m2_0=a0,
                                   min_hgt=h0, max_hgt=h1,
                                   mb_model=mb_mod)

# specify path where to store model diagnostics
diag_path = gdir.get_filepath('model_diagnostics', delete=True)
hef_ref.run_until_and_store(2000, diag_path=diag_path)
ds = utils.compile_run_output([gdir])

ds.to_netcdf(os.path.join(output_dir, 'run_output.nc'))
