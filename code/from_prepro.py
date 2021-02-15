# Locals
from oggm import cfg, utils, workflow
from oggm.core import climate, gis

import oggm_vas as vascaling

# Initialize OGGM and set up the default run parameters
vascaling.initialize()
rgi_version = '62'

# 10 is only for OGGM-VAS, OGGM needs 80 to run
cfg.PARAMS['border'] = 80

# Local working directory (where OGGM will write its output)
WORKING_DIR = utils.gettempdir('moritz')
utils.mkdir(WORKING_DIR, reset=True)
cfg.PATHS['working_dir'] = WORKING_DIR

# RGI glaciers: her you usually pick a region or yours
rgi_ids = ['RGI60-11.00897', 'RGI60-11.00892']
rgi_region = (rgi_ids[0].split('-')[-1]).split('.')[0]

# The important bit: update with the folder with CRU files in it when available
base_url = 'https://cluster.klima.uni-bremen.de/~oggm/gdirs/oggm_v1.4/' \
           'L3-L5_files/RGIV62_fleb_qc3_CRU_pcp2.5'

# Go - get the pre-processed glacier directories
gdirs = workflow.init_glacier_directories(rgi_ids, from_prepro_level=3,
                                          prepro_base_url=base_url,
                                          prepro_rgi_version=rgi_version)

# Start of MY CODE

# define the baseline climate CRU or HISTALP
cfg.PARAMS['baseline_climate'] = 'CRU'
# set the mb hyper parameters accordingly
cfg.PARAMS['prcp_scaling_factor'] = 3
cfg.PARAMS['temp_melt'] = 0
cfg.PARAMS['temp_all_solid'] = 4
cfg.PARAMS['run_mb_calibration'] = False
cfg.PARAMS['use_multiprocessing'] = False
# set minimum ice thickness to include in glacier length computation
# this reduces weird spikes in length records
cfg.PARAMS['min_ice_thick_for_length'] = 0.1

# the bias is defined to be zero during the calibration process,
# which is why we don't use it here to reproduce the results
cfg.PARAMS['use_bias_for_run'] = True

# get and set path to intersect shapefile
intersects_db = utils.get_rgi_intersects_region_file(region=rgi_region)
cfg.set_intersects_db(intersects_db)

# operational run, all glaciers should run
cfg.PARAMS['continue_on_error'] = False

# run vascaling climate tasks
workflow.execute_entity_task(vascaling.local_t_star, gdirs)
# adjust mass balance residual with geodetic observations
vascaling.match_regional_geodetic_mb(gdirs=gdirs, rgi_reg=rgi_region)
# prepare historic "spinup"
workflow.execute_entity_task(vascaling.run_from_climate_data, gdirs=gdirs,
                             ys=2003, ye=2020,
                             output_filesuffix='_historical')
