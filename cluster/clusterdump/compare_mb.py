""" Compare constant mass balance model
---------------------------------------

The VAS model produces highly symmetrical results positive and negative step
changes in climate. This is most likely due to the missing elevation feedback.
In this script I'll compare the constant mass balance models of the scaling
and flowline model.

"""

# import externals libraries
import os
import numpy as np
import pandas as pd
import geopandas as gpd
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


def climate_run_vas(rgi_ids, path=True, temp_biases=[0, +0.5, -0.5],
                    suffixes=['_bias_zero', '_bias_p', '_bias_n'],
                    use_bias_for_run=False, use_default_tstar=True,
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
    kwargs:
        Additional key word arguments for massbalance model.

    Returns
    -------
    Dataset containing yearly values of specific massbalance.

    """
    # assert correct output file suffixes for temp biases
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
    WORKING_DIR = utils.get_temp_dir('_mb_comp')
    OUTPUT_DIR = './'
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
    cfg.PARAMS['prcp_scaling_factor'] = 1.75
    cfg.PARAMS['temp_melt'] = -1.75
    # the bias is defined to be zero during the calibration process,
    # which is why we don't use it here to reproduce the results
    cfg.PARAMS['use_bias_for_run'] = use_bias_for_run

    # operational run, all glaciers should run
    cfg.PARAMS['continue_on_error'] = False

    # read RGI entry for the glaciers as DataFrame
    # containing the outline area as shapefile
    rgidf = utils.get_rgi_glacier_entities(rgi_ids)

    # get and set path to intersect shapefile
    intersects_db = utils.get_rgi_intersects_region_file(region=rgi_region)
    cfg.set_intersects_db(intersects_db)

    # initialize the GlacierDirectory
    gdirs = workflow.init_glacier_directories(rgidf, reset=False, force=True)

    # define the local grid and glacier mask
    workflow.execute_entity_task(gis.define_glacier_region, gdirs)
    workflow.execute_entity_task(gis.glacier_masks, gdirs)
    # process the given climate file
    workflow.execute_entity_task(climate.process_climate_data, gdirs)
    # compute local t* and the corresponding mu*
    if tstar or use_default_tstar:
        # compute mustar from given tstar
        workflow.execute_entity_task(vascaling.local_t_star, gdirs, tstar=tstar, bias=0)
    else:
        # compute mustar from the reference table for the flowline model
        # RGI v6 and HISTALP baseline climate
        ref_df = pd.read_csv(utils.get_demo_file('oggm_ref_tstars_rgi6_histalp.csv'))
        workflow.execute_entity_task(vascaling.local_t_star, gdirs, ref_df=ref_df)

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
        kwargs.setdefault('halfsize', 15)
        kwargs.setdefault('filename', 'climate_historical')
        kwargs.setdefault('input_filesuffix', '')

        ds_ = list()

        for suffix, temp_bias in zip(suffixes, temp_biases):
            # instance mass balance model
            try:
                mb_mod = vascaling.ConstantVASMassBalance(gdir, **kwargs)
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
    kwargs:
        Additional key word arguments for massbalance model.

    Returns
    -------
    Dataset containing yearly values of specific massbalance.
    """

    # assert correct output file suffixes for temp biases
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
    WORKING_DIR = utils.get_temp_dir('_mb_comp')
    OUTPUT_DIR = './'
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
    cfg.PARAMS['prcp_scaling_factor'] = 1.75
    cfg.PARAMS['temp_melt'] = -1.75
    # the bias is defined to be zero during the calibration process,
    # which is why we don't use it here to reproduce the results
    cfg.PARAMS['use_bias_for_run'] = use_bias_for_run

    # operational run, all glaciers should run
    cfg.PARAMS['continue_on_error'] = False

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
        workflow.execute_entity_task(climate.local_t_star, gdirs,
                                     tstar=tstar, bias=0)
    else:
        # compute mustar from the reference table for the flowline model
        # RGI v6 and HISTALP baseline climate
        ref_df = pd.read_csv(utils.get_demo_file('oggm_ref_tstars_rgi6_histalp.csv'))
        workflow.execute_entity_task(climate.local_t_star, gdirs, ref_df=ref_df)
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
        kwargs.setdefault('halfsize', 15)
        kwargs.setdefault('mb_model_class', flowline.ConstantMassBalance)
        kwargs.setdefault('filename', 'climate_historical')
        kwargs.setdefault('input_filesuffix', '')

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


def mb_runs(rgi_ids, tstar=None, path=True):
    """Calls the `climate_run_...` routines for the given RGI IDs,
    combines the resulting datasets for the VAS and flowline model and stores
    it to file.
    """
    vas_ds = climate_run_vas(rgi_ids, tstar=tstar)
    fl_ds = climate_run_fl(rgi_ids, tstar=tstar)

    # concat datasets by evolution balance model
    ds = xr.concat([vas_ds, fl_ds], pd.Index(['vas', 'fl'], name='model'))

    # store to file
    # store to file
    if path:
        if not isinstance(path, str):
            # set default path and filename
            OUTPUT_DIR = os.environ["OUTDIR"]
            path = os.path.join(OUTPUT_DIR, 'mb_output.nc')
        ds.to_netcdf(path)


if __name__ == '__main__':
    rgi_ids = ['RGI60-11.00897']
    nyears = 1
    ds = climate_run_vas(rgi_ids, nyears=nyears)
    climate_run_fl(rgi_ids, nyears=nyears)

