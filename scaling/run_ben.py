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
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt


# import the needed OGGM modules
import oggm
from oggm import cfg
from oggm.utils import get_demo_file
from oggm.tests.funcs import get_test_dir
from oggm.core import gis, climate, centerlines, vascaling

# ---------------------
#  PREPROCESSING TASKS
# ---------------------

# create test directory
testdir = os.path.join(get_test_dir(), 'tmp_ben')
if not os.path.exists(testdir):
    os.makedirs(testdir)
shutil.rmtree(testdir)
os.makedirs(testdir)

# load default parameter file
cfg.initialize()
cfg.set_intersects_db(get_demo_file('rgi_intersect_oetztal.shp'))
cfg.PATHS['working_dir'] = testdir
cfg.PATHS['dem_file'] = get_demo_file('hef_srtm.tif')
# cfg.PATHS['climate_file'] = get_demo_file('histalp_merged_hef.nc')
cfg.PARAMS['border'] = 50
# cfg.PARAMS['run_mb_calibration'] = True
cfg.PARAMS['baseline_climate'] = 'HISTALP'
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
# climate.process_custom_climate_data(gdir)
# climate.process_cru_data(gdir)
climate.process_histalp_data(gdir)

# run center line preprocessing tasks
centerlines.compute_centerlines(gdir)
centerlines.initialize_flowlines(gdir)
centerlines.catchment_area(gdir)
centerlines.catchment_intersections(gdir)
centerlines.catchment_width_geom(gdir)
centerlines.catchment_width_correction(gdir)

# --------------------
#  MASS BALANCE TASKS
# --------------------

# compute local t* and the corresponding mu*
vascaling.local_t_star(gdir)

# instance the mass balance models
mbmod = vascaling.VAScalingMassBalance(gdir)

# ----------------
#  DYNAMICAL PART
# ----------------
# get reference area
a0 = gdir.rgi_area_m2
# get reference year
y0 = gdir.read_pickle('climate_info')['baseline_hydro_yr_0']
# get min and max glacier surface elevation
h0, h1 = vascaling.get_min_max_elevation(gdir)

ben_model = vascaling.VAScalingModel(year_0=y0, area_m2_0=a0,
                                     min_hgt=h0, max_hgt=h1,
                                     mb_model=mbmod)

ben_model.create_start_glacier(area_m2_start=3086069, year_start=1802)

diag_ds = ben_model.run_until_and_store(2003)

store = True
if store:
    # define path and file names
    folder = '/Users/oberrauch/work/master/data/'
    suffix = '_test_start_area'
    names = ['length', 'area', 'volume']
    # combine glacier geometries into DataFrame
    df = diag_ds.to_dataframe()
    df = pd.concat([df['length_m'], df['area_m2'], df['volume_m3']], axis=1)
    df.to_csv(folder+'run'+suffix+'.csv')

plot = True
if plot:
    unit = 'km'
    unit_factor = 1e3
    # define path and file names
    folder = '/Users/oberrauch/work/master/plots/'
    suffix = '_test_start_area'

    plt.figure()
    plt.plot(diag_ds['area_m2']/unit_factor**2)
    plt.title('Hintereis Ferner - Area')
    plt.ylabel('Area [{}$^2$]'.format(unit))
    plt.savefig(folder+'area'+suffix+'.png', bbox_inches='tight')

    plt.figure()
    plt.plot(diag_ds['volume_m3']/unit_factor**3)
    plt.title('Hintereis Ferner - Volume')
    plt.ylabel('Volume [{}$^3$]'.format(unit))
    plt.savefig(folder+'volume'+suffix+'.png', bbox_inches='tight')

    plt.figure()
    plt.plot(diag_ds['length_m']/unit_factor)
    plt.title('Hintereis Ferner - Length')
    plt.ylabel('Length [{}]'.format(unit))
    plt.savefig(folder+'length'+suffix+'.png', bbox_inches='tight')

    plt.show()
