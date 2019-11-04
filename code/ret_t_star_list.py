"""
    Run all necessary preprocessing tasks up to the mu* calibration
    (i.e. `local_mu_star`). Save as... in order to use as template for further
    experiments.

    Date: 07.02.2019
"""
# import externals libs
import os
import shutil
import numpy as np
import pandas as pd
import geopandas as gpd

# import the needed OGGM modules
import oggm
from oggm import cfg
from oggm.utils import get_demo_file, get_rgi_glacier_entities
from oggm.tests.funcs import get_test_dir
from oggm.core import gis, climate, centerlines, vascaling

# ---------------------
#  PREPROCESSING TASKS
# ---------------------

# create test directory
testdir = os.path.join(get_test_dir(), 'tmp_comparison')
if not os.path.exists(testdir):
    os.makedirs(testdir)
shutil.rmtree(testdir)
os.makedirs(testdir)

# load default parameter file
cfg.initialize()
cfg.set_intersects_db(get_demo_file('rgi_intersect_oetztal.shp'))
cfg.PATHS['working_dir'] = testdir
cfg.PATHS['dem_file'] = get_demo_file('hef_srtm.tif')
cfg.PATHS['climate_file'] = get_demo_file('histalp_merged_hef.nc')
cfg.PARAMS['border'] = 10
cfg.PARAMS['run_mb_calibration'] = True
cfg.PARAMS['baseline_climate'] = ''
cfg.PARAMS['use_multiprocessing'] = True

# read the Hintereisferner DEM
hef_file = get_demo_file('Hintereisferner_RGI5.shp')
entity = gpd.read_file(hef_file).iloc[0]

# initialize the GlacierDirectory
gdir = oggm.GlacierDirectory(entity, base_dir=testdir)
# define the local grid and glacier mask
gis.define_glacier_region(gdir, entity=entity)
gis.glacier_masks(gdir)

# process the given climate file
climate.process_custom_climate_data(gdir)

# run center line preprocessing tasks
centerlines.compute_centerlines(gdir)
centerlines.initialize_flowlines(gdir)
centerlines.compute_downstream_line(gdir)
centerlines.compute_downstream_bedshape(gdir)
centerlines.catchment_area(gdir)
centerlines.catchment_intersections(gdir)
centerlines.catchment_width_geom(gdir)
centerlines.catchment_width_correction(gdir)

# get reference mass balance data
mbdf = gdir.get_ref_mb_data()
annual_df = mbdf['ANNUAL BALANCE']
