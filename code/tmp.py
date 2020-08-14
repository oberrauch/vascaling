import pandas as pd
import geopandas as gpd

from oggm import cfg, utils, workflow
from oggm.core import flowline, gis, climate
import oggm_vas as vascaling

rgi_region = '11'
rgi_version = '61'
reset = False

# load default parameter file
cfg.initialize()

# set necessary paths and parameters
wdir = '/Users/oberrauch/work/master/working_directories/histalp'
utils.mkdir(wdir, reset=reset)

cfg.PATHS['working_dir'] = wdir
cfg.PARAMS['rgi_version'] = rgi_version
cfg.PARAMS['baseline_climate'] = 'HISTALP'
cfg.PARAMS['prcp_scaling_factor'] = 1.75
cfg.PARAMS['temp_melt'] = -1.75

cfg.PARAMS['border'] = 10
cfg.PARAMS['use_multiprocessing'] = True
cfg.PARAMS['continue_on_error'] = True

# get RGI dataframe and initialize gdirs
rgi_fpath = utils.get_rgi_region_file(rgi_region, rgi_version, reset=reset)
rgidf = gpd.read_file(rgi_fpath)
# select subregion 1 which is in the HISTALP domain
rgidf = rgidf[rgidf.O2Region == '1']
gdirs = workflow.init_glacier_directories(rgidf)

# # run prepro tasks
# define the local grid and glacier mask
#workflow.execute_entity_task(gis.define_glacier_region, gdirs)
#workflow.execute_entity_task(gis.glacier_masks, gdirs)
# process the given climate file
workflow.execute_entity_task(climate.process_climate_data, gdirs)
