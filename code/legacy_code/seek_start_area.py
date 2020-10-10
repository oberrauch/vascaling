""" First script working on the start area seeking task on the example of
Hintereisferner.

Note: Could be deleted now...
"""

# import externals libs
import os
import shutil
import numpy as np
import geopandas as gpd
from scipy.optimize import minimize_scalar

# import the needed OGGM modules
import oggm
from oggm import cfg
from oggm.utils import get_demo_file, get_rgi_glacier_entities
from oggm.tests.funcs import get_test_dir
from oggm.core import gis, climate, centerlines
from oggm.core import vascaling

# ---------------------
#  PREPROCESSING TASKS
# ---------------------

# create temporary working directory
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
cfg.PARAMS['border'] = 10
# cfg.PARAMS['run_mb_calibration'] = True
cfg.PARAMS['baseline_climate'] = 'HISTALP'
cfg.PARAMS['use_multiprocessing'] = True
# use HistAlp as climate file
# cfg.PATHS['climate_file'] = get_demo_file('histalp_merged_hef.nc')
# set hyper parameters for mass balance calibration
# cfg.PARAMS['baseline_y0'] = 1850
cfg.PARAMS['prcp_scaling_factor'] = 1.75
cfg.PARAMS['temp_melt'] = -1.75

# read the Hintereisferner DEM
hef_file = get_demo_file('Hintereisferner_RGI6.shp')
entity = gpd.read_file(hef_file).iloc[0]

# initialize the GlacierDirectory
gdir = oggm.GlacierDirectory(entity, base_dir=testdir)
# define the local grid and glacier mask
gis.define_glacier_region(gdir, entity=entity)
gis.glacier_masks(gdir)

# process the given climate file
climate.process_histalp_data(gdir)

# run center line preprocessing tasks
centerlines.compute_centerlines(gdir)
centerlines.initialize_flowlines(gdir)
centerlines.compute_downstream_line(gdir)
centerlines.compute_downstream_bedshape(gdir)
centerlines.catchment_area(gdir)
centerlines.catchment_intersections(gdir)
centerlines.catchment_width_geom(gdir)
centerlines.catchment_width_correction(gdir)

# compute the reference t* for the glacier
# given the reference of mass balance measurements
res = vascaling.t_star_from_refmb(gdir)
vascaling.local_t_star(gdir)
t_star, bias = res['t_star'], res['bias']

# --------------------
#  SCALING MODEL
# --------------------

# compute local t* and the corresponding mu*
vascaling.local_t_star(gdir, tstar=t_star, bias=bias)

# instance the mass balance models
ben_mbmod = vascaling.VAScalingMassBalance(gdir)

# get reference area
a0 = gdir.rgi_area_m2
# get reference year
rgi_df = get_rgi_glacier_entities([gdir.rgi_id])
y0 = int(rgi_df.BgnDate.values[0][:4])
# get min and max glacier surface elevation
h0, h1 = vascaling.get_min_max_elevation(gdir)

# set up the glacier model with the values of 2003
hef_ref = vascaling.VAScalingModel(year_0=y0, area_m2_0=a0,
                                   min_hgt=h0, max_hgt=h1,
                                   mb_model=ben_mbmod)


def to_minimize(area_m2_start, model_ref):
    """ Initialize VAS glacier model as copy of the reference model (model_ref)
    and adjust the model to the given starting area (area_m2_start) and
    starting year (1851). Let the model evolve to the same year as the
    reference model. Compute and return the relative absolute area error.

    :param area_m2_start: (float)
    :param model_ref: (vascaling.VAScalingModel)
    :return: (float) relative absolute area error
    """
    # define model
    model_tmp = vascaling.VAScalingModel(year_0=model_ref.year_0,
                                         area_m2_0=model_ref.area_m2_0,
                                         min_hgt=model_ref.min_hgt_0,
                                         max_hgt=model_ref.max_hgt,
                                         mb_model=model_ref.mb_model)
    # scale to desired starting size
    model_tmp.create_start_glacier(area_m2_start, year_start=1851)
    # run and compare, return relative error
    return np.abs(model_tmp.run_and_compare(model_ref))


# define bounds
area_m2_bounds = [0.1, 1e8]
res = minimize_scalar(to_minimize, args=(hef_ref),
                      bounds=area_m2_bounds, method='bounded')
print(res)
