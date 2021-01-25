# Locals
from oggm import cfg, utils, workflow
from oggm.core import climate, gis

import oggm_vas as vascaling

# Initialize OGGM and set up the default run parameters
cfg.initialize()
rgi_version = '62'

# 10 is only for OGGM-VAS, OGGM needs 80 to run
cfg.PARAMS['border'] = 10

# Local working directory (where OGGM will write its output)
WORKING_DIR = utils.gettempdir('moritz')
utils.mkdir(WORKING_DIR, reset=True)
cfg.PATHS['working_dir'] = WORKING_DIR

# RGI glaciers: her you usually pick a region or yours
rgi_ids = ['RGI60-11.00897']
rgi_region = rgi_region = (rgi_ids[0].split('-')[-1]).split('.')[0]

# The important bit: update with the folder with CRU files in it when available
base_url = 'https://cluster.klima.uni-bremen.de/~oggm/gdirs/oggm_v1.4/L2_files/elev_bands/'

# Go - get the pre-processed glacier directories
gdirs = workflow.init_glacier_directories(rgi_ids, from_prepro_level=2,
                                          prepro_base_url=base_url,
                                          prepro_rgi_version=rgi_version)

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
