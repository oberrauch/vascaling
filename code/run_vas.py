""" Hereafter I set up a full run with the volume/area scaling model 
from start to finish. This includes:
a) initialisation and calibration
b) the mass balance model
c) the 'dynamic' model

The results are stored under `run_vas.csv` in the `../data/` directory.
The plots show the temporal evolution of glacier length, surface area
and volume and are store under `length/area/volume_vas.png` in the
`../plots/` directory. 
"""

# import externals libs
import os
import shutil
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt


# import the needed OGGM modules
import oggm
from oggm import cfg, utils
from oggm.core import gis, climate, centerlines, vascaling

# ---------------------
#  PREPROCESSING TASKS
# ---------------------

# define glacier by RGI ID
rgi_id = 'RGI60-11.00737'
# compute RGI region and version from RGI IDs
# assuming all they are all the same
rgi_region = (rgi_id.split('-')[-1]).split('.')[0]
rgi_version = (rgi_id.split('-')[0])[-2:]

# create test directory
wdir = utils.get_temp_dir('tmp_vas')
if not os.path.exists(wdir):
    os.makedirs(wdir)
shutil.rmtree(wdir)
os.makedirs(wdir)

# load default parameter file
cfg.initialize()
# specify some needed parameters and paths
cfg.PATHS['working_dir'] = wdir
cfg.PARAMS['border'] = 50
cfg.PARAMS['baseline_climate'] = 'HISTALP'
cfg.PARAMS['use_multiprocessing'] = True

# get and set path to intersect shapefile
intersects_db = utils.get_rgi_intersects_region_file(region=rgi_region)
cfg.set_intersects_db(intersects_db)

# get RGI entity
entity = utils.get_rgi_glacier_entities([rgi_id])[0]

# initialize the GlacierDirectory
gdir = oggm.GlacierDirectory(entity, base_dir=wdir)
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
tstar = 1927
vascaling.local_t_star(gdir, tstar=tstar)

# instance the mass balance models
mbmod = vascaling.VAScalingMassBalance(gdir)

# ----------------
#  DYNAMICAL PART
# ----------------
# get reference area
a0 = gdir.rgi_area_m2
# get reference year
y0 = gdir.read_pickle('climate_info')['baseline_hydro_yr_0']
y1 = gdir.read_pickle('climate_info')['baseline_hydro_yr_1']
# get min and max glacier surface elevation
h0, h1 = vascaling.get_min_max_elevation(gdir)

vas_model = vascaling.VAScalingModel(year_0=y0, area_m2_0=a0,
                                     min_hgt=h0, max_hgt=h1,
                                     mb_model=mbmod)

# create a 'new' start glacier from a given start area
# vas_model.create_start_glacier(area_m2_start=3086069, year_start=1802)

# run the model
diag_ds = vas_model.run_until_and_store(2003)

# -----------------
#  RESULTS & PLOTS
# -----------------
# store results to file
store = True
if store:
    # define path and file names
    folder = '/Users/oberrauch/work/master/data/'
    suffix = '_vas'
    names = ['length', 'area', 'volume']
    # combine glacier geometries into DataFrame
    df = diag_ds.to_dataframe()
    df = pd.concat([df['length_m'], df['area_m2'], df['volume_m3']], axis=1)
    df.to_csv(folder+'run'+suffix+'.csv')

# create plots and store results to file
plot = True
if plot:
    unit = 'km'
    unit_factor = 1e3
    # define path and file names
    folder = '/Users/oberrauch/work/master/plots/'
    suffix = '_vas'

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
