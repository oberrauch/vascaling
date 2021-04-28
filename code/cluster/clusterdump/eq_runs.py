""" Equilibrium runs
--------------------

This script runs the VAS and flowline model for a single (or more) glacier(s)
with a constant or random massbalance model, performing equilibrium
experiments.

"""

# import internal and externals libraries
import os
import json
import numpy as np
import pandas as pd
import xarray as xr

import logging

log = logging.getLogger('equilibrium-runs')

# import the needed OGGM modules
from oggm import cfg, utils, workflow
from oggm.core import gis, climate, flowline
import oggm_vas as vascaling


def normalize_ds_with_start(ds, store_var_0=False):
    """ Normalize all data variables of the given xarray Dataset
    with their respective first entry. Returns a new xarray.Dataset.

    Parameters
    ----------
    ds: :py:class:`xarray.Dataset`

    Returns
    -------
    Normalized xarray Dataset

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


def equilibrium_run_vas(rgi_ids, use_random_mb=True, path=True,
                        temp_biases=(0, +0.5, -0.5), use_bias_for_run=False,
                        suffixes=('_bias_zero', '_bias_p', '_bias_n'),
                        tstar=None, use_default_tstar=True,
                        vas_c_length_m=None, vas_c_area_m2=None,
                        store_individual_glaciers=True, store_mean_sum=True,
                        prcp_scaling_factor=2.5, temp_melt=-0.5,
                        run_mb_calibration=False, ref_df=None,
                        **kwargs):
    """ The routine runs all steps for the equilibrium experiments using the
    volume/area scaling model:
    - OGGM preprocessing, including initialization, GIS tasks, climate tasks and
      massbalance tasks.
    - Run model for all glaciers with constant (or random) massbalance model
      over 3000 years (default value).
    - Process the model output dataset(s), i.e. normalization, average/sum, ...

    The final dataset containing all results is returned. Given a path it is
    also stored to file.

    Parameters
    ----------
    rgi_ids: array-like
        List of RGI IDs for which the equilibrium experiments are performed.
    use_random_mb: bool, optional, default=True
        Choose between random massbalance model and constant massbalance model.
    path: bool or str, optional, default=True
        If a path is given (or True), the resulting dataset is stored to file.
    temp_biases: array-like, optional, default=(0, +0.5, -0.5)
        List of temperature biases (float, in degC) for the mass balance model.
    suffixes: array-like, optional, default=['_normal', '_bias_p', '_bias_n']
        Descriptive suffixes corresponding to the given temperature biases.
    use_bias_for_run: bool, optional, default=False
        Flag deciding whether or not the mass balance residual is used
    tstar: float, optional, default=None
        'Equilibrium year' used for the mass balance calibration. Using the
        `ref_tstars.csv` table if not supplied.
    use_default_tstar : bool, optional, default=True
        Flag deciding whether or not to use the default ref_tstar.csv list. If
        `False`the `oggm_ref_tstars_rgi6_histalp.csv` reference table is used.
    vas_c_length_m: float, optional, default=None
        Scaling constant for volume/length scaling
    vas_c_area_m2: float, optional, default=None
        Scaling constant for volume/area scaling
    store_mean_sum: bool, optional, default=True
        Flag deciding whether or not to compute and store the average and sum
        over all given glaciers
    store_individual_glaciers : bool, optional, default=True
        Flag deciding whether or not to store the model output for all
        individual glaciers
    prcp_scaling_factor: float, optional, default=2.5
    temp_melt: float, optional, default=-0.5
    run_mb_calibration: float, optional, default=False
    kwargs:
        Additional key word arguments for the `run_random_climate` or
        `run_constant_climate` routines of the vascaling module.

    Returns
    -------
    Dataset containing yearly values of all glacier geometries.

    """
    log.info('Starting VAS equilibrium runs')

    # assert correct output file suffixes for temp biases
    temp_biases = np.atleast_1d(temp_biases)
    suffixes = np.atleast_1d(suffixes)
    if len(temp_biases) != len(suffixes):
        raise RuntimeError("Each given temperature bias must have its "
                           "corresponding suffix")

    # OGGM preprocessing
    # ------------------

    # compute RGI region and version from RGI IDs
    # assuming all they are all the same
    rgi_region = (rgi_ids[0].split('-')[-1]).split('.')[0]
    rgi_version = (rgi_ids[0].split('-')[0])[-2:-1]

    # load default parameter file
    vascaling.initialize()

    # get environmental variables for working and output directories
    WORKING_DIR = os.environ["WORKDIR"]
    OUTPUT_DIR = os.environ["OUTDIR"]
    # create working directory
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
    cfg.PARAMS['prcp_scaling_factor'] = prcp_scaling_factor
    cfg.PARAMS['temp_melt'] = temp_melt
    cfg.PARAMS['run_mb_calibration'] = run_mb_calibration
    # set minimum ice thickness to include in glacier length computation
    # this reduces weird spikes in length records

    # change scaling constants for length and area
    if vas_c_length_m:
        cfg.PARAMS['vas_c_length_m'] = vas_c_length_m
    if vas_c_area_m2:
        cfg.PARAMS['vas_c_area_m2'] = vas_c_area_m2
    # the bias is defined to be zero during the calibration process,
    # which is why we don't use it here to reproduce the results
    cfg.PARAMS['use_bias_for_run'] = use_bias_for_run

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
    if run_mb_calibration and ref_df is not None:
        # compute mustar from given reference table
        log.info('Computing local t* from given reference table')
        workflow.execute_entity_task(vascaling.local_t_star, gdirs,
                                     ref_df=ref_df)
    elif tstar or use_default_tstar:
        # compute mustar from given tstar or reference table
        info = 'Computing local t* from t*' if tstar \
            else 'Computing local t* from default reference table'
        log.info(info)
        workflow.execute_entity_task(vascaling.local_t_star, gdirs,
                                     tstar=tstar, bias=0)
    else:
        # compute mustar from the reference table for the flowline model
        # RGI v6 and HISTALP baseline climate
        log.info('Computing local t* from OGGM RGI6 HISTALP reference table')
        ref_df = pd.read_csv(
            utils.get_demo_file('oggm_ref_tstars_rgi6_histalp.csv'))
        workflow.execute_entity_task(vascaling.local_t_star, gdirs,
                                     ref_df=ref_df)
    # Run model with constant/random mass balance model
    # -------------------------------------------------

    # use t* as center year if not specified otherwise
    kwargs.setdefault('y0', tstar)
    # run for 3000 years if not specified otherwise
    kwargs.setdefault('nyears', 3000)

    if use_random_mb:
        # set random seed to get reproducible results
        kwargs.setdefault('seed', 12)

        # run RandomMassBalance model centered around t*, once without
        # temperature bias and once with positive and negative temperature bias
        # of 0.5 deg C each (per default).
        for suffix, temp_bias in zip(suffixes, temp_biases):
            workflow.execute_entity_task(vascaling.run_random_climate, gdirs,
                                         temperature_bias=temp_bias,
                                         output_filesuffix=suffix, **kwargs)
    else:
        # run ConstantMassBalance model centered around t*, once without
        # temperature bias and once with positive and negative temperature bias
        # of 0.5 deg C each (per default).
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
                                       input_filesuffix=suffix, path=False)
        # add to container
        ds.append(ds_)

    # concat the single output datasets into one,
    # with temperature bias as coordinate
    if ds:
        log.info('Concatenating the output datasets')
        ds = xr.concat(ds, pd.Index(temp_biases, name='temp_bias'))
        # add model type as coordinate
        ds.coords['model'] = 'vas'
        # add mb model type as coordinate
        ds.coords['mb_model'] = 'random' if use_random_mb else 'constant'

        # fill NaN values (which happen for vanished glaciers) with zero
        ds = ds.fillna(0)

        if store_individual_glaciers:
            if store_mean_sum:
                log.info('Computing mean and sum over all glaciers')
                # compute mean and sum over all glaciers
                ds_mean = ds.mean(dim='rgi_id')
                ds_mean.coords['rgi_id'] = 'mean'
                ds_sum = ds.sum(dim='rgi_id')
                ds_sum.coords['rgi_id'] = 'sum'
                # add to datasetv
                ds = xr.concat([ds, ds_mean, ds_sum], dim='rgi_id')
        else:
            log.info('Computing mean and sum over all glaciers')
            # compute mean and sum over all glaciers
            ds_mean = ds.mean(dim='rgi_id')
            ds_mean.coords['rgi_id'] = 'mean'
            ds_sum = ds.sum(dim='rgi_id')
            ds_sum.coords['rgi_id'] = 'sum'
            # add to dataset
            ds = xr.concat([ds_mean, ds_sum], dim='rgi_id')

        log.info('Normalizing glacier geometries with start values')
        # normalize glacier geometries (length/area/volume) with start value
        ds_normal = normalize_ds_with_start(ds)
        # add coordinate to distinguish between normalized and absolute values
        ds.coords['normalized'] = int(False)
        ds_normal.coords['normalized'] = int(True)

        # combine datasets
        ds = xr.concat([ds, ds_normal], 'normalized')

        # store datasets
        if path:
            if path is True:
                mb = 'random' if use_random_mb else 'constant'
                path = os.path.join(OUTPUT_DIR,
                                    'run_output_{}_vas.nc'.format(mb))
            log.info(f'Result are written into {path}')
            ds.to_netcdf(path)

    return ds


def equilibrium_run_fl(rgi_ids, use_random_mb=True, path=True,
                       temp_biases=(0, +0.5, -0.5), use_bias_for_run=False,
                       suffixes=['_bias_zero', '_bias_p', '_bias_n'],
                       store_individual_glaciers=True, store_mean_sum=True,
                       prcp_scaling_factor=1.75, temp_melt=-1.75,
                       run_mb_calibration=False,
                       tstar=None, use_default_tstar=True, **kwargs):
    """ The routine runs all steps for the equilibrium experiments using the
    flowline model. For details see docstring of `sensitivity_run_vas`.

    Parameters
    ----------
    rgi_ids: array-like
        List of RGI IDs for which the equilibrium experiments are performed.
    use_random_mb: bool, optional, default=True
        Choose between random massbalance model and constant massbalance model.
    path: bool or str, optional, default=True
        If a path is given (or True), the resulting dataset is stored to file.
    temp_biases: array-like, optional, default=(0, +0.5, -0.5)
        List of temperature biases (float, in degC) for the mass balance model.
    suffixes: array-like, optional, default=['_normal', '_bias_p', '_bias_n']
        Descriptive suffixes corresponding to the given temperature biases.
    use_bias_for_run: bool, optional, default=False
        Flag deciding whether or not the mass balance residual is used
    tstar: float, optional, default=None
        'Equilibrium year' used for the mass balance calibration. Using the
        `ref_tstars.csv` table if not supplied.
    use_default_tstar : bool, optional, default=True
        Flag deciding whether or not to use the default ref_tstar.csv list. If
        `False`the `oggm_ref_tstars_rgi6_histalp.csv` reference table is used.
    store_mean_sum: bool, optional, default=True
        Flag deciding whether or not to compute and store the average and sum
        over all given glaciers
    store_individual_glaciers : bool, optional, default=True
        Flag deciding whether or not to store the model output for all
        individual glaciers
    prcp_scaling_factor: float, optional, default=1.75
    temp_melt: float, optional, default=-1.75
    run_mb_calibration: float, optional, default=False
    kwargs:
        Additional key word arguments for the `run_random_climate` or
        `run_constant_climate` routines of the vascaling module.

    Returns
    -------
    Dataset containing yearly values of all glacier geometries.

    """
    log.info('Starting flowline equilibrium runs')

    # assert correct output file suffixes for temp biases
    temp_biases = np.atleast_1d(temp_biases)
    suffixes = np.atleast_1d(suffixes)
    if len(temp_biases) != len(suffixes):
        raise RuntimeError("Each given temperature bias must have its "
                           "corresponding suffix")

    # compute RGI region and version from RGI IDs
    # assuming all they are all the same
    rgi_region = (rgi_ids[0].split('-')[-1]).split('.')[0]
    rgi_version = (rgi_ids[0].split('-')[0])[-2:-1]

    # load default parameter file
    cfg.initialize()

    # get environmental variables for working and output directories
    WORKING_DIR = os.environ["WORKDIR"]
    OUTPUT_DIR = os.environ["OUTDIR"]
    # create working directory
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
    cfg.PARAMS['prcp_scaling_factor'] = prcp_scaling_factor
    cfg.PARAMS['temp_melt'] = temp_melt
    cfg.PARAMS['run_mb_calibration'] = run_mb_calibration
    # the bias is defined to be zero during the calibration process,
    # which is why we don't use it here to reproduce the results
    cfg.PARAMS['use_bias_for_run'] = use_bias_for_run
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
    if tstar or use_default_tstar:
        # compute mustar from given tstar or reference table
        info = 'Computing local t* from t*' if tstar \
            else 'Computing local t* from default reference table'
        log.info(info)
        workflow.execute_entity_task(climate.local_t_star, gdirs, tstar=tstar,
                                     bias=0)
    else:
        # compute mustar from the reference table for the flowline model
        # RGI v6 and HISTALP baseline climate
        log.info('Computing local t* from OGGM RGI6 HISTALP reference table')
        ref_df = pd.read_csv(
            utils.get_demo_file('oggm_ref_tstars_rgi6_histalp.csv'))
        workflow.execute_entity_task(climate.local_t_star, gdirs,
                                     ref_df=ref_df)

    workflow.execute_entity_task(climate.mu_star_calibration, gdirs)
    # run inversion tasks
    workflow.inversion_tasks(gdirs)
    # finalize preprocessing
    workflow.execute_entity_task(flowline.init_present_time_glacier, gdirs)

    # use t* as center year if not specified otherwise
    kwargs.setdefault('y0', tstar)
    # run for 3000 years if not specified otherwise
    kwargs.setdefault('nyears', 3000)
    # disregard glaciers exceeding their domain boundaries
    # to not dirsupt the entire run
    kwargs.setdefault('check_for_boundaries', True)

    if use_random_mb:
        # set random seed to get reproducible results
        kwargs.setdefault('seed', 12)

        # run RandomMassBalance model centered around t*, once without
        # temperature bias and once with positive and negative temperature bias
        # of 0.5 deg C each.
        for suffix, temp_bias in zip(suffixes, temp_biases):
            workflow.execute_entity_task(flowline.run_random_climate, gdirs,
                                         temperature_bias=temp_bias,
                                         output_filesuffix=suffix, **kwargs)
    else:
        # run RandomMassBalance model centered around t*, once without
        # temperature bias and once with positive and negative temperature bias
        # of 0.5 deg C each.
        for suffix, temp_bias in zip(suffixes, temp_biases):
            workflow.execute_entity_task(flowline.run_constant_climate, gdirs,
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
                                       input_filesuffix=suffix, path=False)
        # add to container
        ds.append(ds_)

    # concat the single output datasets into one,
    # with temperature bias as coordinate
    if ds:
        log.info('Concatenating the output datasets')
        ds = xr.concat(ds, pd.Index(temp_biases, name='temp_bias'))
        # add model type as coordinate
        ds.coords['model'] = 'fl'
        # add mb model type as coordinate
        ds.coords['mb_model'] = 'random' if use_random_mb else 'constant'

        # fill NaN values (which happen for vanished glaciers) with zero
        ds = ds.fillna(0)

        if store_individual_glaciers:
            if store_mean_sum:
                log.info('Computing mean and sum over all glaciers')
                # compute mean and sum over all glaciers
                ds_mean = ds.mean(dim='rgi_id')
                ds_mean.coords['rgi_id'] = 'mean'
                ds_sum = ds.sum(dim='rgi_id')
                ds_sum.coords['rgi_id'] = 'sum'
                # add to dataset
                ds = xr.concat([ds, ds_mean, ds_sum], dim='rgi_id')
        else:
            log.info('Computing mean and sum over all glaciers')
            # compute mean and sum over all glaciers
            ds_mean = ds.mean(dim='rgi_id')
            ds_mean.coords['rgi_id'] = 'mean'
            ds_sum = ds.sum(dim='rgi_id')
            ds_sum.coords['rgi_id'] = 'sum'
            # add to dataset
            ds = xr.concat([ds_mean, ds_sum], dim='rgi_id')

        log.info('Normalizing glacier geometries with start values')
        # normalize glacier geometries (length/area/volume) with start value
        ds_normal = normalize_ds_with_start(ds)
        # add coordinate to distinguish between normalized and absolute values
        ds.coords['normalized'] = int(False)
        ds_normal.coords['normalized'] = int(True)

        # combine datasets
        ds = xr.concat([ds, ds_normal], 'normalized')

        # store datasets
        if path:
            if path is True:
                mb = 'random' if use_random_mb else 'constant'
                path = os.path.join(OUTPUT_DIR,
                                    'run_output_{}_fl.nc'.format(mb))
            log.info(f'Result are written into {path}')
            ds.to_netcdf(path)

    return ds


def climate_run_vas(rgi_ids, path=True, temp_biases=[0, +0.5, -0.5],
                    suffixes=['_bias_zero', '_bias_p', '_bias_n'],
                    use_bias_for_run=False, use_default_tstar=True,
                    prcp_scaling_factor=2.5, temp_melt=-0.5,
                    run_mb_calibration=False, ref_df=None,
                    tstar=None, nyears=None, **kwargs):
    """Computes 'only' the massbalance in analogy to the `equilibrium_run_...`
    routines, without running the (volume/area scaling) evolution model.

    Dataset containing yearly values of specific mass balance is returned.

    Note: the task is not parallelized, hence it can take long if many glaciers
    are given. TODO: could/should be fixed sometime...
    TODO: add logging information

    Parameters
    ----------
    rgi_ids: array-like
        List of RGI IDs for which the equilibrium experiments are performed.
    path: bool or str, optional, default=True
        If a path is given (or True), the resulting dataset is stored to file.
    temp_biases: array-like, optional, default=(0, +0.5, -0.5)
        List of temperature biases (float, in degC) for the mass balance model.
    suffixes: array-like, optional, default=['_normal', '_bias_p', '_bias_n']
        Descriptive suffixes corresponding to the given temperature biases.
    use_bias_for_run: bool, optional, default=False
        Flag deciding whether or not the mass balance residual is used
    tstar: float, optional, default=None
        'Equilibrium year' used for the mass balance calibration. Using the
        `ref_tstars.csv` table if not supplied.
    use_default_tstar : bool, optional, default=True
        Flag deciding whether or not to use the default ref_tstar.csv list. If
        `False`the `oggm_ref_tstars_rgi6_histalp.csv` reference table is used.
    nyears: int, optional, default=None
        Number of years for which to compute the random mass balance
    prcp_scaling_factor: float, optional, default=1.75
    temp_melt: float, optional, default=-1.75
    run_mb_calibration: float, optional, default=False
    kwargs:
        Additional key word arguments for massbalance model.

    Returns
    -------
    Dataset containing yearly values of specific massbalance.

    """
    # assert correct output file suffixes for temp biases
    temp_biases = np.atleast_1d(temp_biases)
    suffixes = np.atleast_1d(suffixes)
    if len(temp_biases) != len(suffixes):
        raise RuntimeError("Each given temperature bias must have its "
                           "corresponding suffix")

    # compute RGI region and version from RGI IDs
    # assuming all they are all the same
    rgi_region = (rgi_ids[0].split('-')[-1]).split('.')[0]
    rgi_version = (rgi_ids[0].split('-')[0])[-2:-1]

    # load default parameter file
    cfg.initialize()

    # get environmental variables for working and output directories
    WORKING_DIR = os.environ["WORKDIR"]
    OUTPUT_DIR = os.environ["OUTDIR"]
    # create working directory
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
    cfg.PARAMS['prcp_scaling_factor'] = prcp_scaling_factor
    cfg.PARAMS['temp_melt'] = temp_melt
    cfg.PARAMS['run_mb_calibration'] = run_mb_calibration
    # the bias is defined to be zero during the calibration process,
    # which is why we don't use it here to reproduce the results
    cfg.PARAMS['use_bias_for_run'] = use_bias_for_run

    # operational run, all glaciers should run
    cfg.PARAMS['continue_on_error'] = True

    # read RGI entry for the glaciers as DataFrame
    # containing the outline area as shapefile
    rgidf = utils.get_rgi_glacier_entities(rgi_ids)

    # get and set path to intersect shapefile
    intersects_db = utils.get_rgi_intersects_region_file(region=rgi_region)
    cfg.set_intersects_db(intersects_db)

    # initialize the GlacierDirectory
    gdirs = workflow.init_glacier_directories(rgidf, reset=False, force=True)

    # define the local grid and glacier mask
    workflow.execute_entity_task(gis.define_glacier_region, gdirs,
                                 source=None)
    workflow.execute_entity_task(gis.glacier_masks, gdirs)
    # process the given climate file
    workflow.execute_entity_task(climate.process_climate_data, gdirs)
    # compute local t* and the corresponding mu*
    if run_mb_calibration and ref_df is not None:
        # compute mustar from given reference table
        log.info('Computing local t* from given reference table')
        workflow.execute_entity_task(vascaling.local_t_star, gdirs,
                                     ref_df=ref_df)
    elif tstar or use_default_tstar:
        # compute mustar from given tstar or reference table
        info = 'Computing local t* from t*' if tstar \
            else 'Computing local t* from default reference table'
        log.info(info)
        workflow.execute_entity_task(vascaling.local_t_star, gdirs,
                                     tstar=tstar, bias=0)
    else:
        # compute mustar from the reference table for the flowline model
        # RGI v6 and HISTALP baseline climate
        log.info('Computing local t* from OGGM RGI6 HISTALP reference table')
        ref_df = pd.read_csv(
            utils.get_demo_file('oggm_ref_tstars_rgi6_histalp.csv'))
        workflow.execute_entity_task(vascaling.local_t_star, gdirs,
                                     ref_df=ref_df)

    # use t* as center year, even if specified differently
    kwargs['y0'] = tstar
    # run for 10'000 years if not specified otherwise
    if nyears is None:
        nyears = 1e4
    years = np.arange(0, nyears + 1)

    # create dataset
    ds = list()

    # run RandomMassBalance model centered around t*, once without
    # temperature bias and once with positive and negative temperature bias
    # of 0.5 deg C each.
    for gdir in gdirs:
        # set random seed to get reproducible results
        kwargs.setdefault('seed', 12)
        kwargs.setdefault('halfsize', 15)
        kwargs.setdefault('filename', 'climate_historical')
        kwargs.setdefault('input_filesuffix', '')
        kwargs.setdefault('unique_samples', False)

        ds_ = list()

        for suffix, temp_bias in zip(suffixes, temp_biases):
            # instance mass balance model
            try:
                mb_mod = vascaling.RandomVASMassBalance(gdir, **kwargs)
            except:
                # continue with the next glacier or raise exception
                # depending on `continue_on_error` flag
                if cfg.PARAMS['continue_on_error']:
                    continue
                else:
                    raise

            if temp_bias is not None:
                # add given temperature bias to mass balance model
                mb_mod.temp_bias = temp_bias

            # where to store the model output
            diag_path = gdir.get_filepath('model_diagnostics',
                                          filesuffix='_vas',
                                          delete=True)

            # get minimum and maximum glacier elevation
            min_hgt, max_hgt = vascaling.get_min_max_elevation(gdir)
            # create empty container
            spec_mb = list()
            # iterate over all years
            for yr in years:
                spec_mb.append(mb_mod.get_specific_mb(min_hgt, max_hgt, yr))

            # add to dataset
            da = xr.DataArray(spec_mb, dims=('year'), coords={'year': years})
            ds_.append(xr.Dataset({'spec_mb': da}))

        if ds_:
            ds_ = xr.concat(ds_, pd.Index(temp_biases, name='temp_bias'))
            ds_.coords['rgi_id'] = gdir.rgi_id
            ds.append(ds_)

    if ds:
        # combine output from single glaciers into one dataset
        ds = xr.concat(ds, 'rgi_id')

        # store datasets
        if path:
            if path is True:
                path = os.path.join(OUTPUT_DIR, 'mb_output_vas.nc')
            ds.to_netcdf(path)

    # return ds, ds_normal
    return ds


def climate_run_fl(rgi_ids, path=True, temp_biases=[0, +0.5, -0.5],
                   suffixes=['_bias_zero', '_bias_p', '_bias_n'],
                   use_bias_for_run=False, use_default_tstar=True,
                   prcp_scaling_factor=1.75, temp_melt=-1.75,
                   run_mb_calibration=False,
                   tstar=None, nyears=None, **kwargs):
    """Computes 'only' the massbalance in analogy to the `equilibrium_run_...`
    routines, without running the (flowline) evolution model.

    Dataset containing yearly values of specific mass balance is returned.

    Note: the task is not parallelized, hence it can take long if many glaciers
    are given. TODO: could/should be fixed sometime...
    TODO: add logging information

    Parameters
    ----------
    rgi_ids: array-like
        List of RGI IDs for which the equilibrium experiments are performed.
    path: bool or str, optional, default=True
        If a path is given (or True), the resulting dataset is stored to file.
    temp_biases: array-like, optional, default=(0, +0.5, -0.5)
        List of temperature biases (float, in degC) for the mass balance model.
    suffixes: array-like, optional, default=['_normal', '_bias_p', '_bias_n']
        Descriptive suffixes corresponding to the given temperature biases.
    use_bias_for_run: bool, optional, default=False
        Flag deciding whether or not the mass balance residual is used
    tstar: float, optional, default=None
        'Equilibrium year' used for the mass balance calibration. Using the
        `ref_tstars.csv` table if not supplied.
    use_default_tstar : bool, optional, default=True
        Flag deciding whether or not to use the default ref_tstar.csv list. If
        `False`the `oggm_ref_tstars_rgi6_histalp.csv` reference table is used.
    nyears: int, optional, default=None
        Number of years for which to compute the random mass balance
    prcp_scaling_factor: float, optional, default=1.75
    temp_melt: float, optional, default=-1.75
    run_mb_calibration: float, optional, default=False
    kwargs:
        Additional key word arguments for massbalance model.

    Returns
    -------
    Dataset containing yearly values of specific massbalance.
    """

    # assert correct output file suffixes for temp biases
    temp_biases = np.atleast_1d(temp_biases)
    suffixes = np.atleast_1d(suffixes)
    if len(temp_biases) != len(suffixes):
        raise RuntimeError("Each given temperature bias must have its "
                           "corresponding suffix")

    # compute RGI region and version from RGI IDs
    # assuming all they are all the same
    rgi_region = (rgi_ids[0].split('-')[-1]).split('.')[0]
    rgi_version = (rgi_ids[0].split('-')[0])[-2:-1]

    # load default parameter file
    cfg.initialize()

    # get environmental variables for working and output directories
    WORKING_DIR = os.environ["WORKDIR"]
    OUTPUT_DIR = os.environ["OUTDIR"]
    # create working directory
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
    cfg.PARAMS['prcp_scaling_factor'] = prcp_scaling_factor
    cfg.PARAMS['temp_melt'] = temp_melt
    cfg.PARAMS['run_mb_calibration'] = run_mb_calibration
    # the bias is defined to be zero during the calibration process,
    # which is why we don't use it here to reproduce the results
    cfg.PARAMS['use_bias_for_run'] = use_bias_for_run

    # operational run, all glaciers should run
    cfg.PARAMS['continue_on_error'] = True

    # read RGI entry for the glaciers as DataFrame
    # containing the outline area as shapefile
    rgidf = utils.get_rgi_glacier_entities(rgi_ids)

    # get and set path to intersect shapefile
    intersects_db = utils.get_rgi_intersects_region_file(region=rgi_region)
    cfg.set_intersects_db(intersects_db)

    # initialize the GlacierDirectory
    gdirs = workflow.init_glacier_directories(rgidf, reset=False, force=True)

    # run gis tasks
    workflow.gis_prepro_tasks(gdirs)
    # run climate tasks
    workflow.execute_entity_task(climate.process_climate_data, gdirs)
    # compute local t* and the corresponding mu*
    if tstar or use_default_tstar:
        # compute mustar from given tstar
        info = 'Computing local t* from t*' if tstar \
            else 'Computing local t* from default reference table'
        log.info(info)
        workflow.execute_entity_task(climate.local_t_star, gdirs,
                                     tstar=tstar, bias=0)
    else:
        # compute mustar from the reference table for the flowline model
        # RGI v6 and HISTALP baseline climate
        log.info('Computing local t* from OGGM RGI6 HISTALP reference table')
        ref_df = pd.read_csv(
            utils.get_demo_file('oggm_ref_tstars_rgi6_histalp.csv'))
        workflow.execute_entity_task(climate.local_t_star, gdirs,
                                     ref_df=ref_df)
    workflow.execute_entity_task(climate.mu_star_calibration, gdirs)
    # run inversion tasks
    workflow.inversion_tasks(gdirs)
    # finalize preprocessing
    workflow.execute_entity_task(flowline.init_present_time_glacier, gdirs)

    # use t* as center year, even if specified differently
    kwargs['y0'] = tstar
    # run for 3000 years if not specified otherwise
    if nyears is None:
        nyears = 10000
    years = np.arange(0, nyears + 1)

    # create dataset
    ds = list()

    # run RandomMassBalance model centered around t*, once without
    # temperature bias and once with positive and negative temperature bias
    # of 0.5 deg C each.
    for gdir in gdirs:
        # set random seed to get reproducible results
        kwargs.setdefault('seed', 12)
        kwargs.setdefault('halfsize', 15)
        kwargs.setdefault('mb_model_class', flowline.RandomMassBalance)
        kwargs.setdefault('filename', 'climate_historical')
        kwargs.setdefault('input_filesuffix', '')
        kwargs.setdefault('unique_samples', False)

        ds_ = list()

        # open the flowline file if it exists
        try:
            fls = gdir.read_pickle('model_flowlines')
        except:
            # continue with the next glacier or raise exception
            # depending on `continue_on_error` flag
            if cfg.PARAMS['continue_on_error']:
                continue
            else:
                raise

        for suffix, temp_bias in zip(suffixes, temp_biases):
            # instance mass balance model
            try:
                mb_mod = flowline.MultipleFlowlineMassBalance(gdir, **kwargs)
            except:
                # continue with the next glacier or raise exception
                # depending on `continue_on_error` flag
                if cfg.PARAMS['continue_on_error']:
                    continue
                else:
                    raise

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

        if ds_:
            ds_ = xr.concat(ds_, pd.Index(temp_biases, name='temp_bias'))
            ds_.coords['rgi_id'] = gdir.rgi_id
            ds.append(ds_)

    if ds:
        # combine output from single glaciers into one dataset
        ds = xr.concat(ds, 'rgi_id')

        # store dataset to file
        if path:
            if path is True:
                path = os.path.join(OUTPUT_DIR, 'mb_output_fl.nc')
            ds.to_netcdf(path)

    # return ds
    return ds


def eq_runs(rgi_ids, use_random_mb=[True, False], tstar=None, nyears=None,
            use_default_tstar=True, path=True, use_bias_for_run=False,
            store_individual_glaciers=True, store_mean_sum=True,
            run_mb_calibration=False, **kwargs):
    """ Calls the `equilibrium_run_...` routines for the given RGI IDs,
    combines the resulting datasets for the VAS and flowline model and stores
    it to file.
    """
    # perform equilibrium experiments for random and constant climate
    ds = list()
    for r_mb in np.atleast_1d(use_random_mb):
        if not nyears:
            if r_mb:
                # run for 10'000 under random climate
                nyears = 1e4
            else:
                # run foe 1'000 years under constant climate
                nyears = 1e3

        log.info(f"Running a {'random' if r_mb else 'constant'} climate "
                 + f"scenario for {len(rgi_ids):.0f} glaciers "
                 + f"for {nyears:.0f} years")

        vas_ds = equilibrium_run_vas(rgi_ids, use_random_mb=r_mb, tstar=tstar,
                                     use_default_tstar=use_default_tstar,
                                     use_bias_for_run=use_bias_for_run,
                                     path=True, nyears=nyears,
                                     store_individual_glaciers=store_individual_glaciers,
                                     store_mean_sum=store_mean_sum,
                                     run_mb_calibration=run_mb_calibration,
                                     **kwargs)
        fl_ds = equilibrium_run_fl(rgi_ids, use_random_mb=r_mb, tstar=tstar,
                                   use_default_tstar=use_default_tstar,
                                   use_bias_for_run=use_bias_for_run,
                                   run_mb_calibration=run_mb_calibration,
                                   path=True, nyears=nyears,
                                   store_individual_glaciers=store_individual_glaciers,
                                   store_mean_sum=store_mean_sum, **kwargs)
        # concat datasets by model
        ds.append(xr.concat([vas_ds, fl_ds], 'model'))

    # concat datasets by mass balance model
    ds = xr.concat(ds, 'mb_model')

    # store to file
    if path:
        if not isinstance(path, str):
            # set default path and filename
            OUTPUT_DIR = os.environ["OUTDIR"]
            path = os.path.join(OUTPUT_DIR, 'eq_runs.nc')
        ds.to_netcdf(path)


def mb_runs(rgi_ids, tstar=None, use_default_tstar=True, path=True):
    """Calls the `climate_run_...` routines for the given RGI IDs,
    combines the resulting datasets for the VAS and flowline model and stores
    it to file.
    """
    vas_ds = climate_run_vas(rgi_ids, tstar=tstar,
                             use_default_tstar=use_default_tstar)
    fl_ds = climate_run_fl(rgi_ids, tstar=tstar,
                           use_default_tstar=use_default_tstar)

    # concat datasets by evolution balance model
    ds = xr.concat([vas_ds, fl_ds], pd.Index(['vas', 'fl'], name='model'))

    # store to file
    if path:
        if not isinstance(path, str):
            # set default path and filename
            OUTPUT_DIR = os.environ["OUTDIR"]
            path = os.path.join(OUTPUT_DIR, 'mb_output.nc')
        ds.to_netcdf(path)


def vas_eq_runs_ref_df(rgi_ids, ref_df, prcp_scaling_factor, temp_melt,
                       use_random_mb=[True, False], tstar=None, nyears=None,
                       use_default_tstar=True, path=False,
                       use_bias_for_run=False, store_individual_glaciers=True,
                       store_mean_sum=True, **kwargs):
    """ Calls the `equilibrium_run_vas` routines for the given RGI IDs and
    given reference tstar tables, combines the resulting datasets and stores it
    to file.
    """
    # perform equilibrium experiments for random and constant climate
    ds = list()
    for r_mb in np.atleast_1d(use_random_mb):
        if not nyears:
            if r_mb:
                # run for 10'000 under random climate
                nyears = 1e4
            else:
                # run foe 1'000 years under constant climate
                nyears = 1e3

        log.info(f"Running a {'random' if r_mb else 'constant'} climate "
                 + f"scenario for {len(rgi_ids):.0f} glaciers "
                 + f"for {nyears:.0f} years")

        # run vas equilibrium runs
        vas_ds = equilibrium_run_vas(rgi_ids, use_random_mb=r_mb,
                                     tstar=tstar,
                                     use_default_tstar=use_default_tstar,
                                     use_bias_for_run=use_bias_for_run,
                                     path=True, nyears=nyears,
                                     store_individual_glaciers=store_individual_glaciers,
                                     store_mean_sum=store_mean_sum,
                                     prcp_scaling_factor=prcp_scaling_factor,
                                     temp_melt=temp_melt, ref_df=ref_df,
                                     run_mb_calibration=True,
                                     **kwargs)
        # add to container
        ds.append(vas_ds)

    # concat datasets by mass balance model
    ds = xr.concat(ds, 'mb_model')

    # store to file
    if path:
        if not isinstance(path, str):
            # set default path and filename
            OUTPUT_DIR = os.environ["OUTDIR"]
            path = os.path.join(OUTPUT_DIR, 'eq_runs.nc')
        ds.to_netcdf(path)

    return ds


def showcase_glaciers_random_climate():
    """ Run equilibrium experiment for showcase glaciers under random climate
    scenario. """
    # start logger with OGGM settings
    cfg.set_logging_config()

    # get RGI IDs
    fpath = '/home/users/moberrauch/data/showcase_glaciers.csv'
    showcase_glaciers = pd.read_csv(fpath)
    rgi_ids = showcase_glaciers.rgi_id.values

    # start runs
    # mb_runs(rgi_ids)
    eq_runs(rgi_ids, use_random_mb=True, use_default_tstar=False,
            store_mean_sum=False, nyears=2.3e4)


def single_glaciers():
    """ Run equilibrium experiment for showcase glaciers under constant and
    random climate scenario for 1'000 years each."""
    # start logger with OGGM settings
    cfg.set_logging_config()

    # get RGI IDs
    fpath = '/home/users/moberrauch/data/showcase_glaciers.csv'
    showcase_glaciers = pd.read_csv(fpath)
    rgi_ids = showcase_glaciers.rgi_id.values

    # start runs
    eq_runs(rgi_ids, use_default_tstar=False, store_mean_sum=False, nyears=1e3)


def tmp_test():
    """ Run equilibrium experiment for showcase glaciers under constant and
    random climate scenario for 1'000 years each."""
    # start logger with OGGM settings
    cfg.set_logging_config()

    # get RGI IDs
    rgi_ids = ['RGI60-11.00897']

    # start runs
    eq_runs(rgi_ids, use_random_mb=False, use_default_tstar=False, store_mean_sum=False,
            nyears=1e3)


def hef_mb_feedback():
    """ Run equilibrium experiment for Hintereisferner (RGI60-11.00897)
    under constant and random climate scenario for 1'000 years each.
    Thereby the mass balance-elevation feedback of the flowline model is
    ommited using `mb_elev_feedback='never'` as keyword argument."""
    # start logger with OGGM settings
    cfg.set_logging_config()

    # only use Hintereisferner
    rgi_ids = ['RGI60-11.00897']

    # start runs
    eq_runs(rgi_ids, nyears=1e3, use_default_tstar=False, store_mean_sum=False,
            mb_elev_feedback='never')


def histalp_commitment_run():
    """ Run equilibrium experiment for the entire HISTALP domain:
        - constant climate scenario for 1'000 years
        - random climate scenario for 1'000 years
        - each evolution model runs with it's "best fitting" tstar reference table

    """
    # start logger with OGGM settings
    cfg.set_logging_config()

    # get HISTALP RGI IDs
    rgi_ids = pd.read_csv('/home/users/moberrauch/data/histalp_rgi_ids.csv',
                          index_col=0)['RGIId'].values

    # start runs
    eq_runs(rgi_ids, nyears=1e3, use_random_mb=[True, False],
            use_default_tstar=True,
            store_individual_glaciers=False, use_bias_for_run=False)


def histalp_projection_run():
    """ Run equilibrium experiment for the entire HISTALP domain:
        - constant climate scenario for 1'000 years
        - random climate scenario for 1'000 years
        - each evolution model runs with it's "best fitting" tstar reference table

    """
    # start logger with OGGM settings
    cfg.set_logging_config()

    # get HISTALP RGI IDs
    rgi_ids = pd.read_csv('/home/users/moberrauch/data/histalp_rgi_ids.csv',
                          index_col=0)['RGIId'].values

    temp_biases = np.array([0.0, 0.5, 1.0, 1.5, 2.0])
    suffixes = temp_biases.astype(str)

    # start runs
    eq_runs(rgi_ids, nyears=300, use_random_mb=[True, False],
            use_default_tstar=True, temp_biases=temp_biases,
            suffixes=suffixes, y0=1999,
            store_individual_glaciers=False, use_bias_for_run=True)


def compare_vas_ref_dfs():
    """ Run equilibrium experiment for showcase glaciers under constant and
    random climate scenario for 1'000 years each."""
    # start logger with OGGM settings
    cfg.set_logging_config()

    # get RGI IDs
    fpath = '/home/users/moberrauch/data/showcase_glaciers.csv'
    showcase_glaciers = pd.read_csv(fpath)
    rgi_ids = showcase_glaciers.rgi_id.values

    # get reference tables and corresponding
    # mass balance calibration parameters
    ref_dfs = list()
    prcp_scaling_factors = list()
    temp_melt = list()
    dir_path = '/home/users/moberrauch/data/ref_df/'
    for dir in os.listdir(dir_path):
        ref_dfs.append(
            pd.read_csv(os.path.join(dir_path, dir, 'ref_tstars.csv')))
        mb_calib_params = json.load(
            open(os.path.join(dir_path, dir, 'mb_calib_params.json')))
        prcp_scaling_factors.append(mb_calib_params['prcp_scaling_factor'])
        temp_melt.append(mb_calib_params['temp_melt'])

    # create empty container
    vas_ds = list()
    # iterate over all reference tables
    for ref_df, prcpsf, tmelt in zip(ref_dfs, prcp_scaling_factors, temp_melt):
        # equilibrium runs
        vas_ds_ = vas_eq_runs_ref_df(rgi_ids, ref_df, prcpsf, tmelt)
        # define label and add as coordinate
        ref_df_label = f'psf{prcpsf:.2f} tm{tmelt:.2f}'
        vas_ds_.coords['ref_df'] = ref_df_label
        # add to container
        vas_ds.append(vas_ds_)

    # concat into single dataframe
    vas_ds = xr.concat(vas_ds, dim='ref_df')

    # set default path and filename
    OUTPUT_DIR = os.environ["OUTDIR"]
    path = os.path.join(OUTPUT_DIR, 'ref_df_runs.nc')
    # store to file
    vas_ds.to_netcdf(path)


if __name__ == '__main__':
    os.environ[
        'WORKDIR'] = '/Users/oberrauch/work/master/working_directories/test_wd'
    os.environ[
        'OUTDIR'] = '/Users/oberrauch/work/master/working_directories/test_wd'
    tmp_test()
