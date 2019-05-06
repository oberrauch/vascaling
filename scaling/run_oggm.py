"""
    Hereafter I'll try to set up one run with my new model.
    This includes:
    a) initialisation and calibration
    b) the mass balance model
    c) the 'dynamic' model

    Date: 18.01.2019
"""
# import externals libs
import os
import shutil
import geopandas as gpd
import matplotlib.pyplot as plt


# import the needed OGGM modules
import oggm
from oggm.utils import get_demo_file
from oggm.tests.funcs import get_test_dir

# ---------------------
#  PREPROCESSING TASKS
# ---------------------

# create test directory
testdir = os.path.join(get_test_dir(), 'tmp_hef')
if not os.path.exists(testdir):
    os.makedirs(testdir)
shutil.rmtree(testdir)
os.makedirs(testdir)

from oggm import cfg
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

from oggm.core import gis
# initialize the GlacierDirectory
gdir = oggm.GlacierDirectory(entity, base_dir=testdir)
# define the local grid and glacier mask
gis.define_glacier_region(gdir, entity=entity)
gis.glacier_masks(gdir)

from oggm.core import climate
# process the given climate file
climate.process_custom_climate_data(gdir)

from oggm.core import centerlines
# run center line preprocessing tasks
centerlines.compute_centerlines(gdir)
centerlines.initialize_flowlines(gdir)
centerlines.compute_downstream_line(gdir)
centerlines.compute_downstream_bedshape(gdir)
centerlines.catchment_area(gdir)
centerlines.catchment_intersections(gdir)
centerlines.catchment_width_geom(gdir)
centerlines.catchment_width_correction(gdir)

# --------------------
#  MASS BALANCE TASKS
# --------------------

# read reference glacier mass balance data
mbdf = gdir.get_ref_mb_data()
# compute the reference t* for the glacier
# given the reference of mass balance measurements
res = climate.t_star_from_refmb(gdir, mbdf=mbdf['ANNUAL_BALANCE'])
t_star, bias = res['t_star'], res['bias']

# compute local t* and the corresponding mu*
climate.local_t_star(gdir, tstar=t_star, bias=bias)
climate.mu_star_calibration(gdir)

from oggm.core import massbalance
# instance the mass balance models
mb_mod = massbalance.PastMassBalance(gdir)

# -----------
#  INVERSION
# -----------
from oggm.core import inversion
inversion.prepare_for_inversion(gdir)
inversion.mass_conservation_inversion(gdir)
inversion.filter_inversion_output(gdir)

# ----------------
#  DYNAMICAL PART
# ----------------

from oggm.core import flowline
# initialize present time glacier
flowline.init_present_time_glacier(gdir)

# instance flowline model
fls = gdir.read_pickle('model_flowlines')
y0 = gdir.read_pickle('climate_info')['baseline_hydro_yr_0']
fl_mod = flowline.FluxBasedModel(flowlines=fls, mb_model=mb_mod, y0=y0)

# run model and store output as xarray data set
run_ds, diag_ds = fl_mod.run_until_and_store(2003)

# define unit
unit = 'km'
unit_factor = 1e3
# get values
years = diag_ds.hydro_year.values
length = diag_ds.length_m.values / unit_factor
area = diag_ds.area_m2.values / unit_factor**2
volume = diag_ds.volume_m3 / unit_factor**3

store = True
if store:
    # define path and file names
    folder = '/Users/oberrauch/work/master/data/'
    suffix = '_oggm'
    names = ['length', 'area', 'volume']
    # combine glacier geometries into DataFrame
    import numpy as np
    import pandas as pd
    df = pd.DataFrame(np.array([length, area, volume]).T,
                      index=years, columns=[name + suffix for name in names])
    df.to_csv(folder+'run'+suffix+'.csv')

plot = False
if plot:

    # define path and file names
    folder = '/Users/oberrauch/work/master/plots/'
    suffix = '_oggm'

    plt.figure()
    plt.plot(years, area)
    plt.title('Hintereis Ferner - Area')
    plt.ylabel('Area [{}$^2$]'.format(unit))
    plt.savefig(folder+'area'+suffix+'.png', bbox_inches='tight')

    plt.figure()
    plt.plot(years, volume)
    plt.title('Hintereis Ferner - Volume')
    plt.ylabel('Volume [{}$^3$]'.format(unit))
    plt.savefig(folder+'volume'+suffix+'.png', bbox_inches='tight')

    plt.figure()
    plt.plot(years, length)
    plt.title('Hintereis Ferner - Length')
    plt.ylabel('Length [{}]'.format(unit))
    plt.savefig(folder+'length'+suffix+'.png', bbox_inches='tight')
