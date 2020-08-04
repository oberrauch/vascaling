""" The routine runs all steps for the equilibrium experiments using the
volume/area scaling model:
- OGGM preprocessing, including initialization, GIS tasks, climate tasks and
  massbalance tasks.
- Run model for all glaciers with constant (or random) massbalance model
  over 3000 years (default value).
- Process the model output dataset(s), i.e. normalization, average/sum, ...

The final dataset containing all results is returned. Given a path is is
also stored to file.

Parameters
----------
rgi_ids: array-like
    List of RGI IDs for which the equilibrium experiments are performed.
use_random_mb: bool, optional, default=True
    Choose between random massbalance model and constant massbalance model.
use_mean: bool, optional, default=True
    Choose between the mean or summation over all glaciers
path: bool or str, optional, default=True
    If a path is given (or True), the resulting dataset is stored to file.
temp_biases: array-like, optional, default=(0, +0.5, -0.5)
    List of temperature biases (float, in degC) for the mass balance model.
suffixes: array-like, optional, default=['_normal', '_bias_p', '_bias_n']
    Descriptive suffixes corresponding to the given temperature biases.
tstar: float, optional, default=None
    'Equilibrium year' used for the mass balance calibration.
vas_c_length_m: float, optional, default=None
    Scaling constant for volume/length scaling
vas_c_area_m2: float, optional, default=None
    Scaling constant for volume/area scaling
kwargs:
    Additional key word arguments for the `run_random_climate` or
    `run_constant_climate` routines of the vascaling module.

Returns
-------
Dataset containing yearly values of all glacier geometries.

"""


# import externals libraries
import os
import shutil
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

# import the needed OGGM modules
from oggm import cfg, utils, workflow
from oggm.core import gis, climate, flowline, vascaling

rgi_ids = ''
use_random_mb=True
use_mean=True
path=True
temp_biases=(0, +0.5, -0.5)
suffixes=('_normal', '_bias_p', '_bias_n')
tstar=None,
vas_c_length_m=None
vas_c_area_m2=None

kwargs=list()


# assert correct output file suffixes for temp biases
if len(temp_biases) != len(suffixes):
    raise RuntimeError("Each given temperature bias must have its "
                       "corresponding suffix")

# OGGM preprocessing
# ------------------

# compute RGI region and version from RGI IDs
# assuming all they are all the same
rgi_region = (rgi_ids[0].split('-')[-1]).split('.')[0]
rgi_version = (rgi_ids[0].split('-')[0])[-2:]

# load default parameter file
cfg.initialize()

plt.set_xtick(rotation=90)

# create working directory
wdir = '/Users/oberrauch/work/master/working_directories/'
wdir += 'tmp'
if not os.path.exists(wdir):
    os.makedirs(wdir)
# shutil.rmtree(wdir)
# os.makedirs(wdir)
# set path to working directory
cfg.PATHS['working_dir'] = wdir
# set RGI verion and region
cfg.PARAMS['rgi_version'] = rgi_version
# define how many grid points to use around the glacier,
# if you expect the glacier to grow large use a larger border
cfg.PARAMS['border'] = 80
# we use HistAlp climate data
cfg.PARAMS['baseline_climate'] = 'HISTALP'
# set the mb hyper parameters accordingly
cfg.PARAMS['prcp_scaling_factor'] = 1.75
cfg.PARAMS['temp_melt'] = -1.75
# change scaling constants for lenght and area
if vas_c_length_m:
    cfg.PARAMS['vas_c_length_m'] = vas_c_length_m
if vas_c_area_m2:
    cfg.PARAMS['vas_c_area_m2'] = vas_c_area_m2
# the bias is defined to be zero during the calibration process,
# which is why we don't use it here to reproduce the results
cfg.PARAMS['use_bias_for_run'] = False

# read RGI entry for the glaciers as DataFrame
# containing the outline area as shapefile
rgidf = utils.get_rgi_glacier_entities(rgi_ids)

# get and set path to intersect shapefile
intersects_db = utils.get_rgi_intersects_region_file(region=rgi_region)
cfg.set_intersects_db(intersects_db)

# sort by area for more efficient parallel computing
rgidf = rgidf.sort_values('Area', ascending=False)
cfg.PARAMS['use_multiprocessing'] = False

# initialize the GlacierDirectory
gdirs = workflow.init_glacier_regions(rgidf)

# define the local grid and glacier mask
workflow.execute_entity_task(gis.glacier_masks, gdirs)
# process the given climate file
workflow.execute_entity_task(climate.process_histalp_data, gdirs)
# compute local t* and the corresponding mu*
workflow.execute_entity_task(vascaling.local_t_star, gdirs,
                             tstar=tstar, bias=0)

# Run model with constant/random mass balance model
# -------------------------------------------------

# use t* as center year, even if specified differently
kwargs['y0'] = tstar
# run for 3000 years if not specified otherwise
kwargs.setdefault('nyears', 3000)

if use_random_mb:
    # set random seed to get reproducible results
    kwargs.setdefault('seed', 12)

    # run RandomMassBalance model centered around t*, once without
    # temperature bias and once with positive and negative temperature bias
    # of 0.5 °C each (per default).
    for suffix, temp_bias in zip(suffixes, temp_biases):
        workflow.execute_entity_task(vascaling.run_random_climate, gdirs,
                                     temperature_bias=temp_bias,
                                     output_filesuffix=suffix, **kwargs)
else:
    # run ConstantMassBalance model centered around t*, once without
    # temperature bias and once with positive and negative temperature bias
    # of 0.5 °C each (per default).
    for suffix, temp_bias in zip(suffixes, temp_biases):
        workflow.execute_entity_task(vascaling.run_constant_climate, gdirs,
                                     temperature_bias=temp_bias,
                                     output_filesuffix=suffix, **kwargs)

# Process output dataset(s)
# -------------------------

# create empty container
ds = list()
# iterate over all temperature biases/suffixes
for suffix in suffixes:
    # compile the output for each run
    ds_ = utils.compile_run_output(np.atleast_1d(gdirs),
                                   filesuffix=suffix, path=False)
    # add to container
    ds.append(ds_)

# concat the single output datasets into one,
# with temperature bias as coordinate
ds = xr.concat(ds, pd.Index(temp_biases, name='temp_bias'))
# add model type as coordinate
ds.coords['model'] = 'vas'
# add mb model type as coordinate
ds.coords['mb_model'] = 'random' if use_random_mb else 'constant'

# normalize glacier geometries (length/area/volume) with start value
if use_mean:
    # compute average over all glaciers
    ds_normal = normalize_ds_with_start(ds).mean(dim='rgi_id')
    ds = ds.mean(dim='rgi_id')
else:
    # compute sum over all glaciers
    ds_normal = normalize_ds_with_start(ds.sum(dim='rgi_id'))
    ds = ds.sum(dim='rgi_id')

# add coordinate to distinguish between normalized and absolute values
ds.coords['normalized'] = False
ds_normal.coords['normalized'] = True

# combine datasets
ds = xr.concat([ds, ds_normal], 'normalized')

# store datasets
if path:
    if path is True:
        path = list()
        mb = 'random' if use_random_mb else 'constant'
        path.append(os.path.join(cfg.PATHS['working_dir'],
                                 'run_output_{}_vas.nc'.format(mb)))
        # path.append(os.path.join(cfg.PATHS['working_dir'],
        #                          'run_output_{}_vas.nc'.format(mb)))
        # path.append(os.path.join(cfg.PATHS['working_dir'],
        #                          'normalized_output_{}_vas.nc'.format(mb)))
    ds.to_netcdf(path[0])
    # ds_normal.to_netcdf(path[1])


