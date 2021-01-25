# import internal and externals libraries
import os
import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr

import logging

log = logging.getLogger('vas-template')

# import the needed OGGM modules
from oggm import cfg, utils, workflow
from oggm.core import gis, climate, flowline
import oggm_vas as vascaling

log.info('Starting run')

# get all glaciers of RGI region 14 (HIGH MOUNTAIN ASIA)
rgi_region = '14'
rgi_version = '61'
rgidf = gpd.read_file(utils.get_rgi_region_file(region=rgi_region,
                                                 version=rgi_version))

# load default parameter file
vascaling.initialize()

# get path to directories on the CLUSTER - comment/uncomment as necessary
OUTPUT_DIR = os.environ['OUTDIR']
WORKING_DIR = os.environ['WORKDIR']

# set path to working directory
cfg.PATHS['working_dir'] = WORKING_DIR
# set RGI version and region
cfg.PARAMS['rgi_version'] = rgi_version
# define how many grid points to use around the glacier,
# if you expect the glacier to grow large use a larger border
cfg.PARAMS['border'] = 20
# define the baseline cliamte CRU or HISTALP
cfg.PARAMS['baseline_climate'] = 'CRU'
# set the mb hyper parameters accordingly
# cfg.PARAMS['prcp_scaling_factor'] = 2.5
# cfg.PARAMS['temp_melt'] = -0.5
cfg.PARAMS['run_mb_calibration'] = False
# set minimum ice thickness to include in glacier length computation
# this reduces weird spikes in length records
cfg.PARAMS['min_ice_thick_for_length'] = 0.1

# the bias is defined to be zero during the calibration process,
# which is why we don't use it here to reproduce the results
cfg.PARAMS['use_bias_for_run'] = True

# get and set path to intersect shapefile
intersects_db = utils.get_rgi_intersects_region_file(region=rgi_region)
cfg.set_intersects_db(intersects_db)

# sort by area for more efficient parallel computing
rgidf = rgidf.sort_values('Area', ascending=False)
cfg.PARAMS['use_multiprocessing'] = True
# operational run, all glaciers should run
cfg.PARAMS['continue_on_error'] = False

# initialize the GlacierDirectory
gdirs = workflow.init_glacier_directories(rgidf, reset=False, force=True)

# define the local grid and glacier mask
workflow.execute_entity_task(gis.define_glacier_region, gdirs,
                             source=None)
workflow.execute_entity_task(gis.glacier_masks, gdirs)
# process the given climate file
workflow.execute_entity_task(climate.process_climate_data, gdirs)
# compute local t* and the corresponding mu*
workflow.execute_entity_task(vascaling.local_t_star, gdirs)

