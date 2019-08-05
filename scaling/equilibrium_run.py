# Equilibrium runs
# --------------------
# This script runs the VAS model for a single (or more) glacier(s) with a
# constant or random massbalance model, performing equilibrium experiments.
# An additional comparison with the OGGM flowline model is performed.
# TODO

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


def normalize_ds_with_start(ds, store_var_0=False):
    """ Normalize all data variables with their respective first entry.
    Return a new xarray.Dataset.

    Parameters
    ----------
    ds: xarray.Dataset

    Returns
    -------

    """
    # copy dataset
    ds_norm = ds.copy(deep=True)

    # iterate over all data variables
    for var in ds_norm:
        # add information about the initial values
        var_0 = ds_norm[var].isel(time=0)
        if store_var_0:
            ds_norm[var + '_0'] = var_0
        # normalize dataset
        ds_norm[var] /= var_0

    return ds_norm


def equilibrium_run_vas(rgi_ids, use_random_mb=True, use_mean=True, path=True,
                        temp_biases=(0, +0.5, -0.5),
                        suffixes = ['_normal', '_bias_p', '_bias_n'],
                        tstar=None, **kwargs):
    """

    Parameters
    ----------
    rgi_ids
    use_random_mb
    use_mean
    path
    temp_biases
    kwargs

    Returns
    -------

    """
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
    wdir = '/Users/oberrauch/work/master/working_directories/'
    wdir += 'equilibrium_vas_wdir'
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
    # the bias is defined to be zero during the calibration process,
    # which is why we don't use it here to reproduce the results
    cfg.PARAMS['use_bias_for_run'] = False

    # read RGI entry for the glaciers as DataFrame
    # containing the outline area as shapefile
    rgidf = utils.get_rgi_glacier_entities(rgi_ids)

    # get and set path to intersect shapefile
    intersects_db = utils.get_rgi_intersects_region_file(region=rgi_region)
    cfg.set_intersects_db(intersects_db)

    # initialize the GlacierDirectory
    gdirs = workflow.init_glacier_regions(rgidf)

    # define the local grid and glacier mask
    workflow.execute_entity_task(gis.glacier_masks, gdirs)
    # process the given climate file
    workflow.execute_entity_task(climate.process_histalp_data, gdirs)
    # compute local t* and the corresponding mu*
    workflow.execute_entity_task(vascaling.local_t_star, gdirs,
                                 tstar=tstar, bias=0)

    # use t* as center year, even if specified differently
    kwargs['y0'] = tstar
    # run for 3000 years if not specified otherwise
    kwargs.setdefault('nyears', 3000)

    if use_random_mb:
        # set random seed to get reproducible results
        kwargs.setdefault('seed', 12)

        # run RandomMassBalance model centered around t*, once without
        # temperature bias and once with positive and negative temperature bias
        # of 0.5 °C each.
        for suffix, temp_bias in zip(suffixes, temp_biases):
            workflow.execute_entity_task(vascaling.run_random_climate, gdirs,
                                         temperature_bias=temp_bias,
                                         **kwargs, output_filesuffix=suffix)
    else:
        # run ConstantMassBalance model centered around t*, once without
        # temperature bias and once with positive and negative temperature bias
        # of 0.5 °C each.
        for suffix, temp_bias in zip(suffixes, temp_biases):
            workflow.execute_entity_task(vascaling.run_constant_climate, gdirs,
                                         temperature_bias=temp_bias,
                                         **kwargs, output_filesuffix=suffix)

    ds = list()
    for suffix in suffixes:
        # compile the output for each run
        ds_ = utils.compile_run_output(np.atleast_1d(gdirs),
                                       filesuffix=suffix, path=False)
        ds.append(ds_)
    # concat into one dataset with temperature bias as coordinate
    ds = xr.concat(ds, pd.Index(temp_biases, name='temp_bias'))
    # add model type as coordinate
    ds.coords['model'] = 'vas'
    # add mb model type as coordinate
    ds.coords['mb_model'] = 'random' if use_random_mb else 'constant'

    # normalize with start value
    if use_mean:
        ds_normal = normalize_ds_with_start(ds).mean(dim='rgi_id')
        ds = ds.mean(dim='rgi_id')
    else:
        ds_normal = normalize_ds_with_start(ds.sum(dim='rgi_id'))
        ds = ds.sum(dim='rgi_id')

    # add coordinate destiguishing between normalized and absolute values
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

    # return ds, ds_normal
    return ds


def equilibrium_run_fl(rgi_ids, use_random_mb=True, use_mean=True, path=True,
                       temp_biases=(0, +0.5, -0.5),
                       suffixes=['_normal', '_bias_p', '_bias_n'],
                       tstar=None, **kwargs):
    """"""
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
    wdir = '/Users/oberrauch/work/master/working_directories/'
    wdir += 'equilibrium_fl_wdir'
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
    # the bias is defined to be zero during the calibration process,
    # which is why we don't use it here to reproduce the results
    cfg.PARAMS['use_bias_for_run'] = False

    # read RGI entry for the glaciers as DataFrame
    # containing the outline area as shapefile
    rgidf = utils.get_rgi_glacier_entities(rgi_ids)

    # get and set path to intersect shapefile
    intersects_db = utils.get_rgi_intersects_region_file(region=rgi_region)
    cfg.set_intersects_db(intersects_db)

    # initialize the GlacierDirectory
    gdirs = workflow.init_glacier_regions(rgidf)

    # run gis tasks
    workflow.gis_prepro_tasks(gdirs)
    # run climate tasks
    workflow.execute_entity_task(climate.process_histalp_data, gdirs)
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
    kwargs.setdefault('nyears', 3000)

    if use_random_mb:
        # set random seed to get reproducible results
        kwargs.setdefault('seed', 12)

        # run RandomMassBalance model centered around t*, once without
        # temperature bias and once with positive and negative temperature bias
        # of 0.5 °C each.
        for suffix, temp_bias in zip(suffixes, temp_biases):
            workflow.execute_entity_task(flowline.run_random_climate, gdirs,
                                         **kwargs, temperature_bias=temp_bias,
                                         output_filesuffix=suffix)
    else:
        # run RandomMassBalance model centered around t*, once without
        # temperature bias and once with positive and negative temperature bias
        # of 0.5 °C each.
        for suffix, temp_bias in zip(suffixes, temp_biases):
            workflow.execute_entity_task(flowline.run_constant_climate, gdirs,
                                         **kwargs, temperature_bias=temp_bias,
                                         output_filesuffix=suffix)

    ds = list()
    for suffix in suffixes:
        # compile the output for each run
        ds_ = utils.compile_run_output(np.atleast_1d(gdirs),
                                       filesuffix=suffix, path=False)
        ds.append(ds_)
    # concat into one dataset with temperature bias as coordinate
    ds = xr.concat(ds, pd.Index(temp_biases, name='temp_bias'))
    # add model type as coordinate
    ds.coords['model'] = 'fl'
    # add mb model type as coordinate
    ds.coords['mb_model'] = 'random' if use_random_mb else 'constant'

    # normalize with start value
    if use_mean:
        ds_normal = normalize_ds_with_start(ds).mean(dim='rgi_id')
        ds = ds.mean(dim='rgi_id')
    else:
        ds_normal = normalize_ds_with_start(ds.sum(dim='rgi_id'))
        ds = ds.sum(dim='rgi_id')

    # add coordinate destiguishing between normalized and absolute values
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
                                     'run_output_{}_fl.nc'.format(mb)))
            # path.append(os.path.join(cfg.PATHS['working_dir'],
            #                          'run_output_{}_fl.nc'.format(mb)))
            # path.append(os.path.join(cfg.PATHS['working_dir'],
            #                          'normalized_output_{}_fl.nc'.format(mb)))

        ds.to_netcdf(path[0])
        # ds_normal.to_netcdf(path[1])

    # return ds, ds_normal
    return ds


if __name__ == '__main__':
    # define RGI IDs
    rgi_ids = ['RGI60-11.00897']
    # fixate the equilibrium year t*
    tstar = 1927
    # perform equilibrium experiments for random and constant climate
    ds = list()
    for use_random_mb in [True]:
        vas_ds = equilibrium_run_vas(rgi_ids, use_random_mb=use_random_mb,
                                     tstar=tstar, path=False, nyears=1e4)
        fl_ds = equilibrium_run_fl(rgi_ids, use_random_mb=use_random_mb,
                                   tstar=tstar, path=False, nyears=1e4)
        # concat datasets by model
        ds.append(xr.concat([vas_ds, fl_ds], 'model'))

    # concat datasets by mass balance model
    ds = xr.concat(ds, 'mb_model')

    # store to file
    path = '/Users/oberrauch/work/master/data/eq_runs/hef_eq_1e4.nc'
    ds.to_netcdf(path)

