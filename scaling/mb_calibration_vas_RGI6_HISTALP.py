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
from oggm.core import vascaling
from oggm.core import climate
from oggm.core.massbalance import (ConstantMassBalance, PastMassBalance,
                                   MultipleFlowlineMassBalance)


# RGI Version
rgi_version = '61'

# since I'm in the Alps, I'll use HistAlp as baseline climate
baseline = 'HISTALP'

# initialize OGGM and set up the run parameters
cfg.initialize(logging_level='WORKFLOW')

# local paths (where to write the OGGM run output)
dirname = 'OGGM_ref_mb_{}_RGIV{}_OGGM{}'.format(baseline, rgi_version,
                                                oggm.__version__)
WORKING_DIR = utils.gettempdir(dirname, home=True)
utils.mkdir(WORKING_DIR, reset=True)
cfg.PATHS['working_dir'] = WORKING_DIR


# the following code block alters certain parameters from
# the default config file.

# we are running the calibration ourselves
cfg.PARAMS['run_mb_calibration'] = True

# we are using which baseline data?
cfg.PARAMS['baseline_climate'] = baseline

# no need for intersects since this has an effect on the inversion only
cfg.PARAMS['use_intersects'] = False

# use multiprocessing?
cfg.PARAMS['use_multiprocessing'] = True

# set to True for operational runs
cfg.PARAMS['continue_on_error'] = False

if baseline == 'HISTALP':
    # other params: see https://oggm.org/2018/08/10/histalp-parameters/
    cfg.PARAMS['baseline_y0'] = 1850
    cfg.PARAMS['prcp_scaling_factor'] = 1.75
    cfg.PARAMS['temp_melt'] = -1.75


# the next step is to get all the reference glaciers,
# i.e. glaciers with mass balance measurements.


# get the reference glacier ids (they are different for each RGI version)
rgi_dir = utils.get_rgi_dir(version=rgi_version)
df, _ = utils.get_wgms_files()
rids = df['RGI{}0_ID'.format(rgi_version[0])]


# we can't do Antarctica
rids = [rid for rid in rids if not ('-19.' in rid)]

# for HISTALP only RGI reg 11
if baseline == 'HISTALP':
    rids = [rid for rid in rids if '-11.' in rid]

# make a new dataframe with those (this takes a while)
print('Reading the RGI shapefiles...')
rgidf = utils.get_rgi_glacier_entities(rids, version=rgi_version)
print('For RGIV{} we have {} candidate reference '
      'glaciers.'.format(rgi_version, len(rgidf)))

# initialize the glacier regions
gdirs = workflow.init_glacier_regions(rgidf)

# we need to know which period we have data for
print('Process the climate data...')
cfg.PARAMS['continue_on_error'] = True  # Some glaciers are not in Alps
execute_entity_task(tasks.process_histalp_data, gdirs, print_log=False)
cfg.PARAMS['continue_on_error'] = False

# get reference glaciers with mass balance measurements
gdirs = utils.get_ref_mb_glaciers(gdirs)

# keep only these glaciers
rgidf = rgidf.loc[rgidf.RGIId.isin([g.rgi_id for g in gdirs])]

# save to file
rgidf.to_file(os.path.join(WORKING_DIR, 'mb_ref_glaciers.shp'))
print('For RGIV{} and {} we have {} reference glaciers.'.format(rgi_version,
                                                                baseline,
                                                                len(rgidf)))

# sort for more efficient parallel computing
rgidf = rgidf.sort_values('Area', ascending=False)

# initialize glacier directories
gdirs = workflow.init_glacier_regions(rgidf)

# specify needed prepro tasks
task_list = [
    tasks.glacier_masks,
    tasks.compute_centerlines,
    tasks.initialize_flowlines,
    tasks.catchment_area,
    tasks.catchment_intersections,
    tasks.catchment_width_geom,
    tasks.catchment_width_correction,
]
# execute all prepro tasks
for task in task_list:
    execute_entity_task(task, gdirs)

# run climate tasks
vascaling.compute_ref_t_stars(gdirs)
execute_entity_task(vascaling.local_t_star, gdirs)

# we store the associated params
mb_calib = gdirs[0].read_pickle('climate_info')['mb_calib_params']
with open(os.path.join(WORKING_DIR, 'mb_calib_params.json'), 'w') as fp:
    json.dump(mb_calib, fp)
