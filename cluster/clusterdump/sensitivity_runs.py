""" Sensitivity runs
--------------------

This script runs the VAS and flowline model for a single (or more) glacier(s)
with a constant or random massbalance model, performing sensitivity experiments

-> sensitivity to scaling constant
-> sensitivity to time scales

"""

# import externals libraries
import os
import json
import pickle
import numpy as np
import pandas as pd
import xarray as xr

import logging
log = logging.getLogger('sensitivity-runs')

# import the needed OGGM modules
from oggm import cfg, utils, workflow
from oggm.core import gis, climate, flowline
import oggm_vas as vascaling


def recursive_round(li, precision=4):
    """Rounds nested lists and returns in same structure"""
    try:
        return round(li, precision)
    except TypeError:
        return type(li)(recursive_round(x, precision) for x in li)


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


def compute_scaling_params(rgi_ids, path=None):
    """ The routine computes scaling parameters by fitting a linear regression
    to the volume/area and volume/length scatter in log-log space.
    Thereby, the following two cases apply:
    - compute only scaling constants, since scaling exponents have a physical
        basis and should not be changed
    - compute only scaling constants and scaling exponents

    Returns parameters in a 2-level dictionary. The upper level differentiates
    between the two cases, the lower level indicates the parameters.

    Parameters
    ----------
    rgi_ids: array-like
        List of RGI IDs for which the equilibrium experiments are performed.

    Returns
    -------
    Dictionary containing the computed parameters.

    """
    log.info('Starting scaling parameter computation')

    # compute RGI region and version from RGI IDs
    # assuming all they are all the same
    rgi_region = (rgi_ids[0].split('-')[-1]).split('.')[0]
    rgi_version = (rgi_ids[0].split('-')[0])[-2:-1]

    # load default parameter file
    cfg.initialize()

    # get environmental variables for working and output directories
    WORKING_DIR = os.environ["WORKDIR"]
    OUTPUT_DIR = os.environ["OUTDIR"]
    # WORKING_DIR = '/Users/oberrauch/work/master/working_directories/scaling_params/'
    # OUTPUT_DIR = '/Users/oberrauch/work/master/data/scaling_params/'
    # create working directory
    utils.mkdir(WORKING_DIR)
    utils.mkdir(OUTPUT_DIR)
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

    # read RGI entry for the glaciers as DataFrame
    # containing the outline area as shapefile
    rgidf = utils.get_rgi_glacier_entities(rgi_ids)

    # get and set path to intersect shapefile
    intersects_db = utils.get_rgi_intersects_region_file(region=rgi_region)
    cfg.set_intersects_db(intersects_db)
    # the bias is defined to be zero during the calibration process,
    # which is why we don't use it here to reproduce the results
    cfg.PARAMS['use_bias_for_run'] = True

    # sort by area for more efficient parallel computing
    rgidf = rgidf.sort_values('Area', ascending=False)
    cfg.PARAMS['use_multiprocessing'] = True
    # operational run, all glaciers should run
    cfg.PARAMS['continue_on_error'] = False

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

    # create empty dictionary
    params = dict()

    # compute scaling constants for given (fixed) slope
    params['const_only'] = vascaling.get_scaling_constant(gdirs)

    # compute scaling constants and scaling exponent via linear regression
    params['const_expo'] = vascaling.get_scaling_constant_exponent(gdirs)

    # store to file
    if path:
        if not isinstance(path, str):
            # set default path and filename
            path = os.path.join(OUTPUT_DIR, 'scaling_params.json')
        json.dump(params, open(path, 'w'))

    return params


def sensitivity_run_vas(rgi_ids, use_random_mb=False, use_mean=False,
                        path=True, temp_bias=0, tstar=None,
                        use_default_tstar=True, use_bias_for_run=True,
                        scaling_params=[(4.5507, 0.191, 2.2, 1.375)],
                        time_scale_factors=[1],
                        suffixes=['_default'], **kwargs):
    """ The routine runs all steps for the equilibrium experiments using the
    volume/area scaling model (cf. `equilibrium_run_vas`) but for only one
    given temperature bias. However, it is possible to supply a list of
    sensitivity parameters (the scaling constants, and time scale factor) to
    alter the model behavior.
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
    use_mean: bool, optional, default=False
        Choose between the mean or summation over all glaciers
    path: bool or str, optional, default=True
        If a path is given (or True), the resulting dataset is stored to file.
    temp_bias: float, optional, default=0
        Temperature bias (degC) for the mass balance model.
    sensitivity_params: multi-dimensional array-like, optional,
        default=[(4.5507, 0.191, 2.2, 1.375)]
        List containing the scaling constants and scaling exponents for length
        and area scaling as tuples, e.g., (c_l, c_a, q, gamma)
    suffixes: array-like, optional, default=['_default']
        Descriptive suffixes corresponding to the given sensitivity params
    tstar: float, optional, default=None
        'Equilibrium year' used for the mass balance calibration.
    use_default_tstar: bool, optional, default=True
        Flag deciding whether or not to compute mustar from given from reference
        table. Overridden by a given tstar.
    use_bias_for_run: bool, optional, default=True
        Flag deciding whether or not to use the mass balance residual.
    kwargs:
        Additional key word arguments for the `run_random_climate` or
        `run_constant_climate` routines of the vascaling module.

    Returns
    -------
    Dataset containing yearly values of all glacier geometries.

    """
    # assert correct output file suffixes for temp biases
    if len(scaling_params) * len(time_scale_factors) != len(suffixes):
        raise RuntimeError("Each given combination of scaling parameters and "
                           "time scale factor must have its corresponding"
                           "suffix.")

    # OGGM preprocessing
    # ------------------

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
    utils.mkdir(OUTPUT_DIR)
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
    cfg.PARAMS['continue_on_error'] = True

    # initialize the GlacierDirectory
    gdirs = workflow.init_glacier_directories(rgidf, reset=False, force=True)

    # define the local grid and glacier mask
    workflow.execute_entity_task(gis.define_glacier_region, gdirs)
    workflow.execute_entity_task(gis.glacier_masks, gdirs)
    # process the given climate file
    workflow.execute_entity_task(climate.process_climate_data, gdirs)
    # compute local t* and the corresponding mu*
    if tstar or use_default_tstar:
        # compute mustar from given tstar or reference table
        workflow.execute_entity_task(vascaling.local_t_star, gdirs, tstar=tstar, bias=0)
    else:
        # compute mustar from the reference table for the flowline model
        # RGI v6 and HISTALP baseline climate
        ref_df = pd.read_csv(utils.get_demo_file('oggm_ref_tstars_rgi6_histalp.csv'))
        workflow.execute_entity_task(vascaling.local_t_star, gdirs, ref_df=ref_df)

    # Run model with constant/random mass balance model
    # -------------------------------------------------

    # use t* as center year, even if specified differently
    kwargs['y0'] = tstar
    # run for 3000 years if not specified otherwise
    kwargs.setdefault('nyears', 1000)

    # limit parameters to 3 decimal points
    scaling_params = recursive_round(scaling_params, 3)
    time_scale_factors = recursive_round(time_scale_factors, 3)
    # assure that scaling params are handled as tuples (pairs)
    scaling_params_list = np.zeros(len(scaling_params), dtype=object)
    scaling_params_list[:] = scaling_params
    # combine scaling constants, scaling exponents and time scale factor
    # into one iterable array
    sensitivity_params = np.array(np.meshgrid(scaling_params_list,
                                              time_scale_factors)).T
    sensitivity_params = (sensitivity_params
                          .reshape(-1, sensitivity_params.shape[-1]))

    if use_random_mb:
        # set random seed to get reproducible results
        kwargs.setdefault('seed', 12)

        # run RandomMassBalance model centered around t* for each given
        # parameter set
        for suffix, params in zip(suffixes, sensitivity_params):
            cfg.PARAMS['vas_c_length_m'] = params[0][0]
            cfg.PARAMS['vas_c_area_m2'] = params[0][1]
            cfg.PARAMS['vas_q_length'] = params[0][2]
            cfg.PARAMS['vas_gamma_area'] = params[0][3]
            kwargs['time_scale_factor'] = params[1]
            workflow.execute_entity_task(vascaling.run_random_climate, gdirs,
                                         temperature_bias=temp_bias,
                                         output_filesuffix=suffix, **kwargs)
    else:
        # run ConstantMassBalance model centered around t* for each given
        # parameter set
        for suffix, params in zip(suffixes, sensitivity_params):
            cfg.PARAMS['vas_c_length_m'] = params[0][0]
            cfg.PARAMS['vas_c_area_m2'] = params[0][1]
            cfg.PARAMS['vas_q_length'] = params[0][2]
            cfg.PARAMS['vas_gamma_area'] = params[0][3]
            kwargs['time_scale_factor'] = params[1]
            workflow.execute_entity_task(vascaling.run_constant_climate, gdirs,
                                         temperature_bias=temp_bias,
                                         output_filesuffix=suffix, **kwargs)
    # Process output dataset(s)
    # -------------------------

    # create empty container
    ds = list()
    # iterate over all scaling constants
    for i, params in enumerate(scaling_params):
        # create empty container
        ds_ = list()
        # iterate over all time scale factor
        for j, factor in enumerate(time_scale_factors):
            # compile the output for each run
            pos = j + len(time_scale_factors) * i
            ds__ = utils.compile_run_output(np.atleast_1d(gdirs),
                                            filesuffix=suffixes[pos],
                                            path=False)
            # add time scale factor as coordinate
            ds__.coords['time_scale_factor'] = factor
            # add to container
            ds_.append(ds__)

        # concatenate using time scale factor as concat dimension
        ds_ = xr.concat(ds_, dim='time_scale_factor')
        # add scaling constants as coordinate
        params_list = np.zeros(len([params]), dtype=object)
        params_list[:] = [params]
        ds_ = ds_.expand_dims(dim={'scaling_params': params_list})
        # add to container
        ds.append(ds_)

    # concatenate using scaling constants as concat dimension
    ds = xr.concat(ds, dim='scaling_params')

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
        if not isinstance(path, str):
            # set default path and filename
            mb = 'random' if use_random_mb else 'constant'
            path = os.path.join(OUTPUT_DIR, f'run_output_{mb}_vas.nc')
        # write to file
        log.info(f'Writing run output to {path}')
        pickle.dump(ds, open(path, mode='wb'))

    # return ds, ds_normal
    return ds


def histalp_scaling_params():
    """Compute scaling constants and exponents for the HISTALP domain."""
    # start logger with OGGM settings
    cfg.set_logging_config()

    # get HISTALP RGI IDs
    rgi_ids = pd.read_csv('/home/users/moberrauch/data/histalp_rgi_ids.csv', index_col=0)['RGIId'].values
    compute_scaling_params(rgi_ids, path=True)


def histalp_commitment_custom_scaling():
    """Use new found scaling parameters for HISTALP commitment run"""
    # start logger with OGGM settings
    cfg.set_logging_config()

    # get computed scaling parameters
    scaling_params_dict = json.load(open('/home/users/moberrauch/run_output/scaling_params/scaling_params.json'))

    # set scaling parameters
    scaling_params_list = [(4.5507, 0.191, 2.2, 1.375)]  # Global
    const_only = list(scaling_params_dict['const_only'].values())
    const_only.extend([2.2, 1.375])
    scaling_params_list.append(tuple(const_only))  # consts only
    scaling_params_list.append(tuple(list(scaling_params_dict['const_expo'].values())[:4]))  # best fit

    # define file suffixes
    suffixes = ['_default', '_fixed_exp', '_lin_reg']

    # get HISTALP RGI IDs
    rgi_ids = pd.read_csv('/home/users/moberrauch/data/histalp_rgi_ids.csv', index_col=0)['RGIId'].values

    sensitivity_run_vas(rgi_ids=rgi_ids, scaling_params=scaling_params_list,
                        suffixes=suffixes, temp_bias=+0.5)


def histalp_timescale_sensitivity():
    """Use new found scaling parameters for HISTALP commitment run"""
    # start logger with OGGM settings
    cfg.set_logging_config()

    # define file suffixes
    suffixes = ['_default', '_half', '_twice']

    # get HISTALP RGI IDs
    rgi_ids = pd.read_csv('/home/users/moberrauch/data/histalp_rgi_ids.csv', index_col=0)['RGIId'].values

    sensitivity_run_vas(rgi_ids=rgi_ids, time_scale_factors=[1, 0.5, 1.5],
                        suffixes=suffixes, temp_bias=+0.5)


def test():
    """Use new found scaling parameters for HISTALP commitment run"""
    # start logger with OGGM settings
    cfg.set_logging_config()

    # get computed scaling parameters
    scaling_params_dict = json.load(open('/home/users/moberrauch/run_output/scaling_params/scaling_params.json'))

    # set scaling parameters
    scaling_params_list = [tuple(list(scaling_params_dict['const_expo'].values())[:4])]  # Global

    # define file suffixes
    suffixes = ['_lin_reg']

    # get HISTALP RGI IDs
    rgi_ids = pd.read_csv('/home/users/moberrauch/data/histalp_rgi_ids.csv', index_col=0)['RGIId'].values

    sensitivity_run_vas(rgi_ids=rgi_ids, scaling_params=scaling_params_list,
                        suffixes=suffixes, temp_bias=+0.5, nyears=100)
