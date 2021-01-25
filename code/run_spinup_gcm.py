# import internal and externals libraries
import os
import time
import numpy as np
import pandas as pd
import xarray as xr

import matplotlib.pyplot as plt

import logging

log = logging.getLogger('gcm-spinup')

# import the needed OGGM modules
from oggm import cfg, utils, workflow, tasks
from oggm.core import gis, climate, flowline
import oggm_vas as vascaling

# For timing the run
start = time.time()
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
WORKING_DIR = '/Users/oberrauch/work/master/working_directories/gcm-spinup/'
OUTPUT_DIR = '/Users/oberrauch/work/master/data/gcm-spinup/'
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
# cfg.PARAMS['prcp_scaling_factor'] = 2.5
# cfg.PARAMS['temp_melt'] = -0.5
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

# -----------------
# PREPRO DONE
# -----------------

# Additional climate file (CESM)
cfg.PATHS['cesm_temp_file'] = utils.get_demo_file('cesm.TREFHT.160001-200512'
                                                  '.selection.nc')
cfg.PATHS['cesm_precc_file'] = utils.get_demo_file('cesm.PRECC.160001-200512'
                                                   '.selection.nc')
cfg.PATHS['cesm_precl_file'] = utils.get_demo_file('cesm.PRECL.160001-200512'
                                                   '.selection.nc')
workflow.execute_entity_task(tasks.process_cesm_data, gdirs)

# Run the last 200 years with the default starting point (current glacier)
# and CESM data as input
workflow.execute_entity_task(vascaling.run_from_climate_data, gdirs,
                             climate_filename='gcm_data',
                             ys=1801, ye=2000,
                             output_filesuffix='_no_spinup')

# Run the spinup simulation - t* climate with a cold temperature bias
workflow.execute_entity_task(vascaling.run_constant_climate, gdirs,
                             nyears=100, bias=0, temperature_bias=-0.5,
                             output_filesuffix='_spinup')
# Run a past climate run based on this spinup
workflow.execute_entity_task(vascaling.run_from_climate_data, gdirs,
                             climate_filename='gcm_data',
                             ys=1801, ye=2000,
                             init_model_filesuffix='_spinup',
                             output_filesuffix='_with_spinup')

# Compile output
log.info('Compiling output')
utils.compile_glacier_statistics(gdirs)
ds1 = utils.compile_run_output(gdirs, input_filesuffix='_no_spinup')
ds2 = utils.compile_run_output(gdirs, input_filesuffix='_with_spinup')

# Log
m, s = divmod(time.time() - start, 60)
h, m = divmod(m, 60)
log.info('OGGM is done! Time needed: %d:%02d:%02d' % (h, m, s))

# Plot
f, ax = plt.subplots(figsize=(9, 4))
(ds1.volume.sum(dim='rgi_id') * 1e-9).plot(ax=ax, label='No spinup')
(ds2.volume.sum(dim='rgi_id') * 1e-9).plot(ax=ax, label='With spinup')
ax.set_ylabel('Volume (km$^3$)')
ax.set_xlabel('')
ax.set_title('Hintereisferner volume under CESM forcing')
plt.legend()
plt.tight_layout()
plt.show()
