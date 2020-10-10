""" Compares the results of the original Matlab model to the results of the
OGGM VAS model, during the iterative process of seeking the start area.

The *.mat files are produced by the original Matlab code from Marzeion, et. al.
They contain yearly values for glacier geometries (length, area and volume),
time scales for length and area change (tau_L and tau_A). One file corresponds
to one iteration step.

The corresponding values of the OGGM VAS model are computed during the
iteration. Figures showing length, area and volume, and time scales for length
and area change as functions of time for both models are plotted and displayed.

"""
# import externals libs
import os
import shutil
from scipy.io import loadmat
import matplotlib.pyplot as plt

# import the needed OGGM modules
import oggm
from oggm import cfg, utils
from oggm.tests.funcs import get_test_dir
from oggm.core import gis, climate, centerlines, vascaling

# create temporary working directory
testdir = os.path.join(get_test_dir(), 'start_area')
if not os.path.exists(testdir):
    os.makedirs(testdir)
shutil.rmtree(testdir)
os.makedirs(testdir)

# load default parameter file
cfg.initialize()
# set path to working directory
cfg.PATHS['working_dir'] = testdir
# load and set path to intersects
path = utils.get_rgi_intersects_region_file('11', version='6')
cfg.set_intersects_db(path)

# change some default parameters
cfg.PARAMS['border'] = 50
cfg.PARAMS['baseline_climate'] = 'CRU'
cfg.PARAMS['use_multiprocessing'] = True

# get RGI entity for Hintereisferner
rgi_id = 'RGI60-11.00897'
entity = utils.get_rgi_glacier_entities([rgi_id]).iloc[0]

# initialize the GlacierDirectory
gdir = oggm.GlacierDirectory(entity, base_dir=testdir)

# define the local grid and glacier mask
gis.define_glacier_region(gdir, entity=entity)
gis.glacier_masks(gdir)

# process the given climate file
climate.process_cru_data(gdir)

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
# get reference area (from RGI entry)
a0 = gdir.rgi_area_m2
# get reference year (start of climate records)
y0 = gdir.read_pickle('climate_info')['baseline_hydro_yr_0']
# get min and max glacier surface elevation (based on RGI outline)
h0, h1 = vascaling.get_min_max_elevation(gdir)

# initialize iteration counter variable
k = 1

# start iteration process of "finding the start area"
while True:
    # specify path to *.mat file and see if file exists
    mat_file = '../data/start_area_results/{:s}_iteration_{:02d}.mat'.format(rgi_id, k)
    if not os.path.isfile(mat_file):
        # break loop if no file found
        break

    # read mat file
    mat_file = loadmat(mat_file)
    # increment counter
    k += 1

    # initialize VAS model
    vas_model = vascaling.VAScalingModel(year_0=y0, area_m2_0=a0,
                                         min_hgt=h0, max_hgt=h1,
                                         mb_model=mbmod)

    # get start area from *.mat file
    area_m2_start = mat_file['A_pre'][0, 0] * 1e6

    # init model new start area
    vas_model.create_start_glacier(area_m2_start, year_start=1902)
    diag_ds = vas_model.run_until_and_store(2003)

    # prepare parameters for plotting
    mat_years = mat_file['years'][0]
    vas_years = diag_ds.hydro_year
    mat_params = ['A_pre', 'L_pre', 'V_pre', 'mb_modeled_pre', 'tau_L', 'tau_A']
    vas_params = ['area_m2', 'length_m', 'volume_m3', 'spec_mb', 'tau_l', 'tau_a']
    labels = ['Surface area [km$^2$]', 'Glacier length [km]',
              'Glacier volume [km3]', 'Specific mass balance [mm w.e. yr$^{-1}$]',
              'Length change response time [yr]', 'Area change response time [yr]']
    factors = [1e6, 1e3, 1e9, 1, 1, 1]
    # plot all parameters
    for mat, vas, label, f in zip(mat_params, vas_params, labels, factors):
        plt.figure()
        plt.plot(mat_years, mat_file[mat][0, :102], label='Marzeion Matlab')
        plt.plot(vas_years, diag_ds[vas] / f, label='OGGM VAS')
        plt.legend()
        plt.xlabel('')
        plt.ylabel(label)
        plt.title('Hintereisferner with CRU data')

    plt.show()

    # wait for keyboard input
    input('Compare next iteration step...')
