# import externals libs
import os
import shutil
import numpy as np
import pandas as pd
import geopandas as gpd
from scipy.optimize import minimize, minimize_scalar

# import the needed OGGM modules
import oggm
from oggm import cfg, utils
from oggm.utils import get_rgi_glacier_entities
from oggm.tests.funcs import get_test_dir
from oggm.core import gis, climate, centerlines
from oggm.core import vascaling

# load parametere file
cfg.initialize()

# get/downlaod the rgi entity including the outline shapefile
rgidf = get_rgi_glacier_entities(['RGI60-11.01270'])
rgi_entity = rgidf.iloc[0]
rgidf.plot()

# specify the working directory and define the glacier directory
cfg.PATHS['working_dir'] = './'
gdir = oggm.GlacierDirectory(rgidf.loc[1269])

# get the path to the DEM file (will download if necessary)
dem = utils.get_topo_file(gdir.cenlon, gdir.cenlat)
print('DEM source: {}, path to DEM file: {}'.format(dem[1], dem[0][0]))
# set path in config file
cfg.PATHS['dem_file'] = dem[0][0]
cfg.PARAMS['border'] = 10
cfg.PARAMS['use_intersects'] = False

# gis tasks
gis.define_glacier_region(gdir, entity=rgi_entity)
gis.glacier_masks(gdir)

# climate tasks
cfg.PARAMS['baseline_climate'] = 'HISTALP'
climate.process_histalp_data(gdir)
