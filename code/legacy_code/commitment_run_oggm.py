""" Run the OGGM flowline model for all glaciers of the Rofental Basin
under a random climate based on the HistAlp climate between 1984 and 
2014. The random model is initialized with the same seed (12), the get
reproducible results.

The run output is stored with the suffix `_commitment` in the working
directory (`../working_directories/commitment_run_oggm/`) as netCDF
file.
"""

# import external libs
import logging
import os
import geopandas as gpd
import shapely.geometry as shpg

# import oggm modules
import oggm.cfg as cfg
from oggm import utils, workflow, tasks
from oggm.core import flowline

# For timing the run
import time
start = time.time()

# Module logger
log = logging.getLogger(__name__)

# Initialize OGGM and set up the default run parameters
cfg.initialize(logging_level='WORKFLOW')
rgi_version = '61'
rgi_region = '11'

# Local working directory (where OGGM will write its output)
wdir = '../working_directories/commitment_run_oggm/'
utils.mkdir(wdir, reset=True)
cfg.PATHS['working_dir'] = wdir

# We use intersects
path = utils.get_rgi_intersects_region_file(rgi_region,
    version=rgi_version)
cfg.set_intersects_db(path)

# RGI file
path = utils.get_rgi_region_file(rgi_region, version=rgi_version)
rgidf = gpd.read_file(path)

# Get the Rofental Basin file
path = utils.get_demo_file('rofental_hydrosheds.shp')
basin = gpd.read_file(path)

# Take all glaciers in the Rofental Basin
in_bas = [basin.geometry.contains(shpg.Point(x, y))[0] for
          (x, y) in zip(rgidf.CenLon, rgidf.CenLat)]
rgidf = rgidf.loc[in_bas]

# Sort for more efficient parallel computing
rgidf = rgidf.sort_values('Area', ascending=False)
cfg.PARAMS['use_multiprocessing'] = True

log.workflow('Starting OGGM run')
log.workflow('Number of glaciers: {}'.format(len(rgidf)))

# Go - initialize glacier directories
gdirs = workflow.init_glacier_regions(rgidf)

# run prepro tasks
cfg.PARAMS['baseline_climate'] = 'HISTALP'
cfg.PARAMS['prcp_scaling_factor'] = 1.75
cfg.PARAMS['temp_melt'] = -1.75

# execute GIS tasks
workflow.gis_prepro_tasks(gdirs)
# execute climate tasks
workflow.climate_tasks(gdirs)
# execute inversion tasks
workflow.inversion_tasks(gdirs)
# convert inversion flowlines into model flowlines
workflow.execute_entity_task(flowline.init_present_time_glacier, gdirs)

# Random climate representative for the recent climate (1984-2014)
# This is a kind of "commitment" run
nyears = 500
workflow.execute_entity_task(tasks.run_random_climate, gdirs,
                             nyears=nyears, y0=1999, seed=12,
                             output_filesuffix='_commitment')

# Write the compiled output
utils.compile_run_output(gdirs, filesuffix='_commitment')

# Log
m, s = divmod(time.time() - start, 60)
h, m = divmod(m, 60)
log.workflow('OGGM is done! Time needed: %d:%02d:%02d' % (h, m, s))
