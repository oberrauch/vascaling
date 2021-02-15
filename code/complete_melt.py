# import internal and externals libraries
import os
import numpy as np
import pandas as pd
import xarray as xr

import logging

log = logging.getLogger('vas-template')

# import the needed OGGM modules
from oggm import cfg, utils, workflow
from oggm.core import gis, climate, flowline
import oggm_vas as vascaling

log.info('Starting run')

# specify glaciers by RGI IDs (INPUT)
rgi_ids = ['RGI60-11.00897']

# compute RGI region and version from RGI IDs
# assuming all they are all the same
rgi_region = (rgi_ids[0].split('-')[-1]).split('.')[0]
rgi_version = (rgi_ids[0].split('-')[0])[-2:-1]

# load default parameter file
vascaling.initialize()

# get LOCAL environmental variables for working and output directories
WORKING_DIR = '/Users/oberrauch/work/master/working_directories/complete_melt/'
OUTPUT_DIR = '/Users/oberrauch/work/master/data/complete_melt/'
# create working directory
utils.mkdir(WORKING_DIR)
utils.mkdir(OUTPUT_DIR)

# get path to directories on the CLUSTER - comment/uncomment as necessary
# OUTPUT_DIR = os.environ['OUTDIR']
# WORKING_DIR = os.environ['WORKDIR']


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
cfg.PARAMS['temp_melt'] = 0
cfg.PARAMS['prcp_scaling_factor'] = 3
cfg.PARAMS['temp_all_solid'] = 4
cfg.PARAMS['run_mb_calibration'] = False
# set minimum ice thickness to include in glacier length computation
# this reduces weird spikes in length records
cfg.PARAMS['min_ice_thick_for_length'] = 0.1

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
cfg.PARAMS['use_multiprocessing'] = False
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

vascaling.run_constant_climate(gdirs[0], nyears=300, temperature_bias=15,
                               output_filesuffix='_complete_melt')


