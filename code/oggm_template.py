# import internal and externals libraries
import os
import numpy as np
import pandas as pd
import xarray as xr

import logging

log = logging.getLogger('oggm-template')

# import the needed OGGM modules
from oggm import cfg, utils, workflow
from oggm.core import gis, climate, flowline

log.info('Starting run')

# specify glaciers by RGI IDs (INPUT)
rgi_ids = ['RGI60-11.00897']

# compute RGI region and version from RGI IDs
# assuming all RGI IDs are from within one version and region
rgi_region = (rgi_ids[0].split('-')[-1]).split('.')[0]
rgi_version = (rgi_ids[0].split('-')[0])[-2:-1]

# load default parameter file
cfg.initialize()

# specify path to working and output directories (INPUT)
WORKING_DIR = '/Users/oberrauch/work/master/working_directories/oggm_run/'
OUTPUT_DIR = '/Users/oberrauch/work/master/data/oggm_run/'

# create working directory
utils.mkdir(WORKING_DIR)
utils.mkdir(OUTPUT_DIR)
# set path to working directory
cfg.PATHS['working_dir'] = WORKING_DIR
# set RGI version and region
cfg.PARAMS['rgi_version'] = rgi_version
# define how many grid points to use around the glacier,
# if you expect the glacier to grow large use a larger border
cfg.PARAMS['border'] = 10
# use default climate (=CRU)
# for HISTALP uncomment the following lines (INPUT)
# cfg.PARAMS['baseline_climate'] = 'HISTALP'
# cfg.PARAMS['prcp_scaling_factor'] = 1.75
# cfg.PARAMS['temp_melt'] = -1.75
# cfg.PARAMS['run_mb_calibration'] = False

# the bias is defined to be zero during the calibration process,
# which is why we don't use it here to reproduce the results
cfg.PARAMS['use_bias_for_run'] = True
# set minimum ice thickness to include in glacier length computation
# this reduces weird spikes in length records
cfg.PARAMS['min_ice_thick_for_length'] = 0.1

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
cfg.PARAMS['continue_on_error'] = True

# initialize the GlacierDirectory
gdirs = workflow.init_glacier_directories(rgidf, reset=False, force=True)

# run gis tasks
workflow.gis_prepro_tasks(gdirs)
# run climate tasks
workflow.execute_entity_task(climate.process_climate_data, gdirs)
# compute local t* and the corresponding mu*
workflow.execute_entity_task(climate.local_t_star, gdirs)
workflow.execute_entity_task(climate.mu_star_calibration, gdirs)
# run inversion tasks
workflow.inversion_tasks(gdirs)
# finalize preprocessing
workflow.execute_entity_task(flowline.init_present_time_glacier, gdirs)

# # use t* as center year if not specified otherwise
# kwargs = list()
# kwargs.setdefault('y0', tstar=None)
# # run for 3000 years if not specified otherwise
# kwargs.setdefault('nyears', 3000)
# # disregard glaciers exceeding their domain boundaries
# # to not dirsupt the entire run
# kwargs.setdefault('check_for_boundaries', True)

# if use_random_mb:
#     # set random seed to get reproducible results
#     kwargs.setdefault('seed', 12)
#
#     # run RandomMassBalance model centered around t*, once without
#     # temperature bias and once with positive and negative temperature bias
#     # of 0.5 deg C each.
#     for suffix, temp_bias in zip(suffixes, temp_biases):
#         workflow.execute_entity_task(flowline.run_random_climate, gdirs,
#                                      temperature_bias=temp_bias,
#                                      output_filesuffix=suffix, **kwargs)
# else:
#     # run RandomMassBalance model centered around t*, once without
#     # temperature bias and once with positive and negative temperature bias
#     # of 0.5 deg C each.
#     for suffix, temp_bias in zip(suffixes, temp_biases):
#         workflow.execute_entity_task(flowline.run_constant_climate, gdirs,
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
#     ds.coords['model'] = 'fl'
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
#             # add to dataset
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
#                                 'run_output_{}_fl.nc'.format(mb))
#         log.info(f'Result are written into {path}')
#         ds.to_netcdf(path)
#
