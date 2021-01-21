# import internal and externals libraries
import os
import numpy as np
import pandas as pd
import xarray as xr

import logging

log = logging.getLogger('vas-template')

# import the needed OGGM modules
from oggm import cfg, utils, workflow
from oggm.core import gis, climate, flowline
import oggm_vas as vascaling

log.info('Starting run')

# specify glaciers by RGI IDs (INPUT)
rgi_ids = ['RGI60-11.00897']

# compute RGI region and version from RGI IDs
# assuming all they are all the same
rgi_region = (rgi_ids[0].split('-')[-1]).split('.')[0]
rgi_version = (rgi_ids[0].split('-')[0])[-2:-1]

# load default parameter file
vascaling.initialize()

# get environmental variables for working and output directories
WORKING_DIR = '/Users/oberrauch/work/master/working_directories/vas_run/'
OUTPUT_DIR = '/Users/oberrauch/work/master/data/vas_run/'

# create working directory
utils.mkdir(WORKING_DIR)
utils.mkdir(OUTPUT_DIR)
# set path to working directory
cfg.PATHS['working_dir'] = WORKING_DIR
# set RGI version and region
cfg.PARAMS['rgi_version'] = rgi_version
# define how many grid points to use around the glacier,
# if you expect the glacier to grow large use a larger border
cfg.PARAMS['border'] = 20
# define the baseline cliamte CRU or HISTALP
cfg.PARAMS['baseline_climate'] = 'CRU'
# set the mb hyper parameters accordingly
# cfg.PARAMS['prcp_scaling_factor'] = 2.5
# cfg.PARAMS['temp_melt'] = -0.5
cfg.PARAMS['run_mb_calibration'] = False
# set minimum ice thickness to include in glacier length computation
# this reduces weird spikes in length records
cfg.PARAMS['min_ice_thick_for_length'] = 0.1

# the bias is defined to be zero during the calibration process,
# which is why we don't use it here to reproduce the results
cfg.PARAMS['use_bias_for_run'] = True

# read RGI entry for the glaciers as DataFrame
# containing the outline area as shapefile
rgidf = utils.get_rgi_glacier_entities(rgi_ids)

# get and set path to intersect shapefile
intersects_db = utils.get_rgi_intersects_region_file(region=rgi_region)
cfg.set_intersects_db(intersects_db)

# sort by area for more efficient parallel computing
rgidf = rgidf.sort_values('Area', ascending=False)
cfg.PARAMS['use_multiprocessing'] = True
# operational run, all glaciers should run
cfg.PARAMS['continue_on_error'] = False

# initialize the GlacierDirectory
gdirs = workflow.init_glacier_directories(rgidf, reset=False, force=True)

# define the local grid and glacier mask
workflow.execute_entity_task(gis.define_glacier_region, gdirs,
                             source=None)
workflow.execute_entity_task(gis.glacier_masks, gdirs)
# process the given climate file
workflow.execute_entity_task(climate.process_climate_data, gdirs)
# compute local t* and the corresponding mu*
workflow.execute_entity_task(vascaling.local_t_star, gdirs)

# elif tstar or use_default_tstar:
#     # compute mustar from given tstar or reference table
#     info = 'Computing local t* from t*' if tstar \
#         else 'Computing local t* from default reference table'
#     log.info(info)
#     workflow.execute_entity_task(vascaling.local_t_star, gdirs,
#                                  tstar=tstar, bias=0)
# else:
#     # compute mustar from the reference table for the flowline model
#     # RGI v6 and HISTALP baseline climate
#     log.info('Computing local t* from OGGM RGI6 HISTALP reference table')
#     ref_df = pd.read_csv(
#         utils.get_demo_file('oggm_ref_tstars_rgi6_histalp.csv'))
#     workflow.execute_entity_task(vascaling.local_t_star, gdirs,
#                                  ref_df=ref_df)
# # Run model with constant/random mass balance model
# # -------------------------------------------------
#
# # use t* as center year if not specified otherwise
# kwargs.setdefault('y0', tstar)
# # run for 3000 years if not specified otherwise
# kwargs.setdefault('nyears', 3000)
#
# if use_random_mb:
#     # set random seed to get reproducible results
#     kwargs.setdefault('seed', 12)
#
#     # run RandomMassBalance model centered around t*, once without
#     # temperature bias and once with positive and negative temperature bias
#     # of 0.5 deg C each (per default).
#     for suffix, temp_bias in zip(suffixes, temp_biases):
#         workflow.execute_entity_task(vascaling.run_random_climate, gdirs,
#                                      temperature_bias=temp_bias,
#                                      output_filesuffix=suffix, **kwargs)
# else:
#     # run ConstantMassBalance model centered around t*, once without
#     # temperature bias and once with positive and negative temperature bias
#     # of 0.5 deg C each (per default).
#     for suffix, temp_bias in zip(suffixes, temp_biases):
#         workflow.execute_entity_task(vascaling.run_constant_climate, gdirs,
#                                      temperature_bias=temp_bias,
#                                      output_filesuffix=suffix, **kwargs)
#
# # Process output dataset(s)
# # -------------------------
#
# # create empty container
# ds = list()
# # iterate over all temperature biases/suffixes
# for suffix in suffixes:
#     # compile the output for each run
#     ds_ = utils.compile_run_output(np.atleast_1d(gdirs),
#                                    input_filesuffix=suffix, path=False)
#     # add to container
#     ds.append(ds_)
#
# # concat the single output datasets into one,
# # with temperature bias as coordinate
# if ds:
#     log.info('Concatenating the output datasets')
#     ds = xr.concat(ds, pd.Index(temp_biases, name='temp_bias'))
#     # add model type as coordinate
#     ds.coords['model'] = 'vas'
#     # add mb model type as coordinate
#     ds.coords['mb_model'] = 'random' if use_random_mb else 'constant'
#
#     # fill NaN values (which happen for vanished glaciers) with zero
#     ds = ds.fillna(0)
#
#     if store_individual_glaciers:
#         if store_mean_sum:
#             log.info('Computing mean and sum over all glaciers')
#             # compute mean and sum over all glaciers
#             ds_mean = ds.mean(dim='rgi_id')
#             ds_mean.coords['rgi_id'] = 'mean'
#             ds_sum = ds.sum(dim='rgi_id')
#             ds_sum.coords['rgi_id'] = 'sum'
#             # add to datasetv
#             ds = xr.concat([ds, ds_mean, ds_sum], dim='rgi_id')
#     else:
#         log.info('Computing mean and sum over all glaciers')
#         # compute mean and sum over all glaciers
#         ds_mean = ds.mean(dim='rgi_id')
#         ds_mean.coords['rgi_id'] = 'mean'
#         ds_sum = ds.sum(dim='rgi_id')
#         ds_sum.coords['rgi_id'] = 'sum'
#         # add to dataset
#         ds = xr.concat([ds_mean, ds_sum], dim='rgi_id')
#
#     log.info('Normalizing glacier geometries with start values')
#     # normalize glacier geometries (length/area/volume) with start value
#     ds_normal = normalize_ds_with_start(ds)
#     # add coordinate to distinguish between normalized and absolute values
#     ds.coords['normalized'] = int(False)
#     ds_normal.coords['normalized'] = int(True)
#
#     # combine datasets
#     ds = xr.concat([ds, ds_normal], 'normalized')
#
#     # store datasets
#     if path:
#         if path is True:
#             mb = 'random' if use_random_mb else 'constant'
#             path = os.path.join(OUTPUT_DIR,
#                                 'run_output_{}_vas.nc'.format(mb))
#         log.info(f'Result are written into {path}')
#         ds.to_netcdf(path)
#
# return ds