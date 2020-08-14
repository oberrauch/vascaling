""" Equilibrium runs
--------------------

This script runs the VAS and flowline model for a single (or more) glacier(s)
with a constant or random massbalance model, performing equilibrium
experiments.

"""

# import externals libraries
import os
import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr

# import the needed OGGM modules
from oggm import cfg, utils, workflow
from oggm.core import gis, climate, flowline
import oggm_vas as vascaling

# set method arguments here
path = True
temp_biases = [0, +0.5, -0.5],
suffixes = ['_bias_zero', '_bias_p', '_bias_n'],
tstar = None
nyears = None
kwargs = dict()

# define RGI IDS
rgi_ids = ['RGI60-11.01519', 'RGI60-11.01511', 'RGI60-11.01480']


# assert correct output file suffixes for temp biases
if len(temp_biases) != len(suffixes):
    raise RuntimeError("Each given temperature bias must have its "
                       "corresponding suffix")

# compute RGI region and version from RGI IDs
# assuming all they are all the same
rgi_region = (rgi_ids[0].split('-')[-1]).split('.')[0]
rgi_version = (rgi_ids[0].split('-')[0])[-2:]

# load default parameter file
cfg.initialize()

# create working directory
# WORKING_DIR = os.environ["WORKDIR"]
WORKING_DIR = utils.get_temp_dir('minimal_working_example', home=True)
utils.mkdir(WORKING_DIR)
# set path to working directory
cfg.PATHS['working_dir'] = WORKING_DIR
# set RGI version and region
cfg.PARAMS['rgi_version'] = rgi_version
# define how many grid points to use around the glacier,
# if you expect the glacier to grow large use a larger border
cfg.PARAMS['border'] = 120
# we use HistAlp climate data
cfg.PARAMS['baseline_climate'] = 'HISTALP'
# set the mb hyper parameters accordingly
cfg.PARAMS['prcp_scaling_factor'] = 1.75
cfg.PARAMS['temp_melt'] = -1.75
# the bias is defined to be zero during the calibration process,
# which is why we don't use it here to reproduce the results
cfg.PARAMS['use_bias_for_run'] = False

# operational run, all glaciers should run
cfg.PARAMS['continue_on_error'] = True

# read RGI entry for the glaciers as DataFrame
# containing the outline area as shapefile
rgidf = utils.get_rgi_glacier_entities(rgi_ids)

# get and set path to intersect shapefile
intersects_db = utils.get_rgi_intersects_region_file(region=rgi_region)
cfg.set_intersects_db(intersects_db)

# initialize the GlacierDirectory
gdirs = workflow.init_glacier_directories(rgidf, reset=True, force=True)

# run gis tasks
workflow.gis_prepro_tasks(gdirs)
# run climate tasks
workflow.execute_entity_task(climate.process_climate_data, gdirs)
workflow.execute_entity_task(climate.local_t_star, gdirs,
                             tstar=tstar, bias=0)
workflow.execute_entity_task(climate.mu_star_calibration, gdirs)
# run inversion tasks
workflow.inversion_tasks(gdirs)
# finalize preprocessing
workflow.execute_entity_task(flowline.init_present_time_glacier, gdirs)

# use t* as center year, even if specified differently
kwargs['y0'] = tstar
# run for 3000 years if not specified otherwise
if nyears is None:
    nyears = 1000
years = np.arange(0, nyears + 1)

# create dataset
ds = list()

# run RandomMassBalance model centered around t*, once without
# temperature bias and once with positive and negative temperature bias
# of 0.5 °C each.
for gdir in gdirs:
    # set random seed to get reproducible results
    kwargs.setdefault('seed', 12)
    kwargs.setdefault('halfsize', 15)
    kwargs.setdefault('mb_model_class', flowline.RandomMassBalance)
    kwargs.setdefault('filename', 'climate_historical')
    kwargs.setdefault('input_filesuffix', '')
    kwargs.setdefault('unique_samples', False)

    ds_ = list()
    try:
        fls = gdir.read_pickle('model_flowlines')
    except FileNotFoundError:
        if cfg.PARAMS['continue_on_error']:
            continue
        else:
            raise

    for suffix, temp_bias in zip(suffixes, temp_biases):
        # instance mass balance model
        mb_mod = flowline.MultipleFlowlineMassBalance(gdir, **kwargs)

        if temp_bias is not None:
            # add given temperature bias to mass balance model
            mb_mod.temp_bias = temp_bias

        # create empty container
        spec_mb = list()
        # iterate over all years
        for yr in years:
            spec_mb.append(mb_mod.get_specific_mb(fls=fls, year=yr))

        # add to dataset
        da = xr.DataArray(spec_mb, dims=('year'), coords={'year': years})
        ds_.append(xr.Dataset({'spec_mb': da}))

    ds_ = xr.concat(ds_, pd.Index(temp_biases, name='temp_bias'))
    ds_.coords['rgi_id'] = gdir.rgi_id
    ds.append(ds_)

if ds:
    ds = xr.concat(ds, 'rgi_id')

# store datasets
if path and ds:
    if path is True:
        path = os.path.join(cfg.PATHS['working_dir'], 'mb_output_fl.nc')
    ds.to_netcdf(path)
