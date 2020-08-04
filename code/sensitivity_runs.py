""" Equilibrium runs
--------------------

This script runs the VAS model for a single (or more) glacier(s)
with a constant or random massbalance model, performing sensitivity
experiments (varying the scaling constants and time scales)

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


def sensitivity_run_vas_old(rgi_ids, use_random_mb=False, use_mean=True,
                            path=True, temp_bias=0, tstar=None,
                            sensitivity_params=[[(4.5507, 0.191), 1]], suffixes=[''],
                            **kwargs):
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
    use_mean: bool, optional, default=True
        Choose between the mean or summation over all glaciers
    path: bool or str, optional, default=True
        If a path is given (or True), the resulting dataset is stored to file.
    temp_bias: float, optional, default=0
        Temperature bias (degC) for the mass balance model.
    sensitivity_params: multi-dimensional array-like, optional,
        default=[[(4.5507, 0.191), 1]]
        list containing the parameters which are to be varied in the following
        order: float tuple with length and area scaling constant, float as time
        scale factor
    suffixes: array-like, optional, default=['']
        Descriptive suffixes corresponding to the given sensitivity params
    tstar: float, optional, default=None
        'Equilibrium year' used for the mass balance calibration.
    kwargs:
        Additional key word arguments for the `run_random_climate` or
        `run_constant_climate` routines of the vascaling module.

    Returns
    -------
    Dataset containing yearly values of all glacier geometries.

    """
    # assert correct output file suffixes for temp biases
    if len(sensitivity_params) != len(suffixes):
        raise RuntimeError("Each given parameter set must have its "
                           "corresponding suffix")

    # OGGM preprocessing
    # ------------------

    # compute RGI region and version from RGI IDs
    # assuming all they are all the same
    rgi_region = (rgi_ids[0].split('-')[-1]).split('.')[0]
    rgi_version = (rgi_ids[0].split('-')[0])[-2:]

    # load default parameter file
    cfg.initialize()

    # create working directory
    wdir = '/Users/oberrauch/work/master/working_directories/'
    wdir += 'sensitivity_vas_wdir'
    if not os.path.exists(wdir):
        os.makedirs(wdir)
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

    # Run model with constant/random mass balance model
    # -------------------------------------------------

    # use t* as center year, even if specified differently
    kwargs['y0'] = tstar
    # run for 3000 years if not specified otherwise
    kwargs.setdefault('nyears', 3000)

    if use_random_mb:
        # set random seed to get reproducible results
        kwargs.setdefault('seed', 12)

        # run RandomMassBalance model centered around t* for each given
        # parameter set
        for suffix, params in zip(suffixes, sensitivity_params):
            cfg.PARAMS['vas_c_length_m'] = params[0]
            cfg.PARAMS['vas_c_area_m2'] = params[1]
            kwargs['time_scale_factor'] = params[2]
            workflow.execute_entity_task(vascaling.run_random_climate, gdirs,
                                         temperature_bias=temp_bias,
                                         output_filesuffix=suffix, **kwargs)
    else:
        # run ConstantMassBalance model centered around t* for each given
        # parameter set
        for suffix, params in zip(suffixes, sensitivity_params):
            cfg.PARAMS['vas_c_length_m'] = params[0][0]
            cfg.PARAMS['vas_c_area_m2'] = params[0][1]
            kwargs['time_scale_factor'] = params[1]
            workflow.execute_entity_task(vascaling.run_constant_climate, gdirs,
                                         temperature_bias=temp_bias,
                                         output_filesuffix=suffix, **kwargs)
    # Process output dataset(s)
    # -------------------------

    # create empty container
    ds = list()
    # iterate over all temperature biases/suffixes
    for suffix, params in zip(suffixes, sensitivity_params):
        # compile the output for each run
        ds_ = utils.compile_run_output(np.atleast_1d(gdirs),
                                       filesuffix=suffix, path=False)
        # add sensitivity parameters as coordinates
        ds_.coords['length_scaling_const'] = params[0][0]
        ds_.coords['area_scaling_const'] = params[0][1]
        ds_.coords['time_scale_factor'] = params[1]
        # add to container
        ds.append(ds_)

    # concat the single output datasets into one, using 'sensitivity_params'
    # as name fot the new concatenate dimension
    ds = xr.combine_nested(ds, )
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
        ds.to_netcdf(path)
        # ds_normal.to_netcdf(path[1])

    # return ds, ds_normal
    return ds


def sensitivity_run_vas(rgi_ids, use_random_mb=False, use_mean=True,
                        path=True, temp_bias=0, tstar=None,
                        scaling_constants=[(4.5507, 0.191)],
                        time_scale_factors=[1],
                        suffixes=[''], **kwargs):
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
    use_mean: bool, optional, default=True
        Choose between the mean or summation over all glaciers
    path: bool or str, optional, default=True
        If a path is given (or True), the resulting dataset is stored to file.
    temp_bias: float, optional, default=0
        Temperature bias (degC) for the mass balance model.
    sensitivity_params: multi-dimensional array-like, optional,
        default=[[(4.5507, 0.191), 1]]
        list containing the parameters which are to be varied in the following
        order: float tuple with length and area scaling constant, float as time
        scale factor
    suffixes: array-like, optional, default=['']
        Descriptive suffixes corresponding to the given sensitivity params
    tstar: float, optional, default=None
        'Equilibrium year' used for the mass balance calibration.
    kwargs:
        Additional key word arguments for the `run_random_climate` or
        `run_constant_climate` routines of the vascaling module.

    Returns
    -------
    Dataset containing yearly values of all glacier geometries.

    """
    # assert correct output file suffixes for temp biases
    if len(scaling_constants) * len(time_scale_factors) != len(suffixes):
        raise RuntimeError("Each given combination of scaling constants and "
                           "time scale factor must have its corresponding"
                           "suffix.")

    # OGGM preprocessing
    # ------------------

    # compute RGI region and version from RGI IDs
    # assuming all they are all the same
    rgi_region = (rgi_ids[0].split('-')[-1]).split('.')[0]
    rgi_version = (rgi_ids[0].split('-')[0])[-2:]

    # load default parameter file
    cfg.initialize()

    # create working directory
    wdir = '/Users/oberrauch/work/master/working_directories/'
    wdir += 'sensitivity_vas_wdir'
    if not os.path.exists(wdir):
        os.makedirs(wdir)
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

    # Run model with constant/random mass balance model
    # -------------------------------------------------

    # use t* as center year, even if specified differently
    kwargs['y0'] = tstar
    # run for 3000 years if not specified otherwise
    kwargs.setdefault('nyears', 3000)

    def recursive_round(li, precision=3):
        """Rounds nested lists and returns in same structure"""
        try:
            return round(li, precision)
        except TypeError:
            return type(li)(recursive_round(x, precision) for x in li)

    # limit parameters to 3 decimal points
    scaling_constants = recursive_round(scaling_constants, 3)
    time_scale_factors = recursive_round(time_scale_factors, 3)
    # assure that scaling constants are handled as tuples (pairs)
    scaling_constants_list = np.zeros(len(scaling_constants), dtype=object)
    scaling_constants_list[:] = scaling_constants
    # combine scaling constants and time scale factor into one iterable array
    sensitivity_params = np.array(np.meshgrid(scaling_constants_list,
                                              time_scale_factors)).T
    sensitivity_params = (sensitivity_params
                          .reshape(-1, sensitivity_params.shape[-1]))

    if use_random_mb:
        # set random seed to get reproducible results
        kwargs.setdefault('seed', 12)

        # run RandomMassBalance model centered around t* for each given
        # parameter set
        # for i, consts in enumerate(scaling_constants):
        #     cfg.PARAMS['vas_c_length_m'] = consts[0]
        #     cfg.PARAMS['vas_c_area_m2'] = consts[1]
        #     for j, factor in enumerate(time_scale_factors):
        #         kwargs['time_scale_factor'] = factor
        #         workflow.execute_entity_task(vascaling.run_random_climate,
        #                                      gdirs, temperature_bias=temp_bias,
        #                                      output_filesuffix=suffixes[2*i+j],
        #                                      **kwargs)

        for suffix, params in zip(suffixes, sensitivity_params):
            cfg.PARAMS['vas_c_length_m'] = params[0]
            cfg.PARAMS['vas_c_area_m2'] = params[1]
            kwargs['time_scale_factor'] = params[2]
            workflow.execute_entity_task(vascaling.run_random_climate, gdirs,
                                         temperature_bias=temp_bias,
                                         output_filesuffix=suffix, **kwargs)
    else:
        # run ConstantMassBalance model centered around t* for each given
        # parameter set
        for suffix, params in zip(suffixes, sensitivity_params):
            cfg.PARAMS['vas_c_length_m'] = params[0][0]
            cfg.PARAMS['vas_c_area_m2'] = params[0][1]
            kwargs['time_scale_factor'] = params[1]
            workflow.execute_entity_task(vascaling.run_constant_climate, gdirs,
                                         temperature_bias=temp_bias,
                                         output_filesuffix=suffix, **kwargs)
    # Process output dataset(s)
    # -------------------------

    # create empty container
    ds = list()
    # iterate over all scaling constants
    for i, consts in enumerate(scaling_constants):
        # create empty container
        ds_ = list()
        # iterate over all time scale factor
        for j, factor in enumerate(time_scale_factors):
            # compile the output for each run
            pos = j+len(time_scale_factors)*i
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
        consts_list = np.zeros(len([consts]), dtype=object)
        consts_list[:] = [consts]
        ds_ = ds_.expand_dims(dim={'scaling_const': consts_list})
        # ds_.coords['area_scaling_const'] = consts[1]
        # add to container
        ds.append(ds_)

    # concatenate using scaling constants as concat dimension
    # ds = xr.concat(ds, dim='area_scaling_const')
    ds = xr.concat(ds, dim='scaling_const')

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
            mb = 'random' if use_random_mb else 'constant'
            path = os.path.join(cfg.PATHS['working_dir'],
                                     'run_output_{}_vas.nc'.format(mb))
            # path.append(os.path.join(cfg.PATHS['working_dir'],
            #                          'run_output_{}_vas.nc'.format(mb)))
            # path.append(os.path.join(cfg.PATHS['working_dir'],
            #                          'normalized_output_{}_vas.nc'.format(mb)))
        import pickle
        pickle.dump(ds, open(path, mode='wb'))
        # ds_normal.to_netcdf(path[1])

    # return ds, ds_normal
    return ds


def time_scale_sensitivity():
    """

    Returns
    -------

    """
    # define RGI IDs
    rgi_ids = ['RGI60-11.00897']
    # fixate the equilibrium year t*
    tstar = 1927

    # define scaling constants
    # scaling_const_list = [(1.5038582086557708, 0.24399290770672957)]  # HEF
    scaling_const_list = [(4.5507, 0.191)]  # Global
    scaling_const = np.zeros(len(scaling_const_list), dtype=object)
    scaling_const[:] = scaling_const_list
    # define time scale factor
    time_scale_factor = np.unique(np.concatenate([1/np.arange(1, 11),
                                                  np.arange(1, 11)]))
    sensitivity_params = np.array(np.meshgrid(scaling_const,
                                              time_scale_factor)).T
    sensitivity_params = (sensitivity_params
                          .reshape(-1, sensitivity_params.shape[-1]))
    # specify file suffixes
    suffixes = ['_{:2f}_{:2f}_{:.2f}'.format(p[0][0], p[0][1], p[1])
                for p in sensitivity_params]

    # specify file storage location
    path = '/Users/oberrauch/work/master/data/eq_runs/hef_time_sensitivity.nc'

    # perform equilibrium experiments for constant climate
    sensitivity_run_vas(rgi_ids, use_random_mb=False, tstar=tstar,
                        temp_bias=-0.5, sensitivity_params=sensitivity_params,
                        suffixes=suffixes, path=path, nyears=3e4)


def scaling_const_sensitivity():
    """

    Returns
    -------

    """
    # define RGI IDs
    rgi_ids = ['RGI60-11.00897']
    # fixate the equilibrium year t*
    tstar = 1927

    # define scaling constants
    scaling_const_list = list((1.5038582086557708, 0.24399290770672957))  # HEF
    scaling_const_list.append((4.5507, 0.191))  # Global
    scaling_const = np.zeros(len(scaling_const_list), dtype=object)
    scaling_const[:] = scaling_const_list
    # define time scale factor
    time_scale_factor = [1]
    sensitivity_params = np.array(np.meshgrid(scaling_const,
                                              time_scale_factor)).T
    sensitivity_params = (sensitivity_params
                          .reshape(-1, sensitivity_params.shape[-1]))
    # specify file suffixes
    suffixes = ['_{:2f}_{:2f}_{:.2f}'.format(p[0][0], p[0][1], p[1])
                for p in sensitivity_params]

    # specify file storage location
    path = '/Users/oberrauch/work/master/data/eq_runs/hef_const_sensitivity.nc'

    # perform equilibrium experiments for constant climate
    sensitivity_run_vas(rgi_ids, use_random_mb=False, tstar=tstar,
                        temp_bias=-0.5, sensitivity_params=sensitivity_params,
                        suffixes=suffixes, path=path, nyears=7e3)


def sensitivity():
    """

    Returns
    -------

    """
    # define RGI IDs
    rgi_ids = ['RGI60-11.00897']
    # fixate the equilibrium year t*
    tstar = 1927

    # define scaling constants
    scaling_const_list = [(1.5038582086557708, 0.24399290770672957)]  # HEF
    scaling_const_list.append((4.5507, 0.191))  # Global
    scaling_const_list.append((5, 0.15))  # Test
    scaling_constants = np.zeros(len(scaling_const_list), dtype=object)
    scaling_constants[:] = scaling_const_list
    # define time scale factor
    time_scale_factors = [0.5, 1, 2]
    sensitivity_params = np.array(np.meshgrid(scaling_constants,
                                              time_scale_factors)).T
    sensitivity_params = (sensitivity_params
                          .reshape(-1, sensitivity_params.shape[-1]))
    # specify file suffixes
    suffixes = ['_{:.3f}_{:.3f}_{:.3f}'.format(p[0][0], p[0][1], p[1])
                for p in sensitivity_params]

    # specify file storage location
    path = '/Users/oberrauch/work/master/data/eq_runs/hef_sensitivity_growth.nc'

    # perform equilibrium experiments for constant climate
    sensitivity_run_vas(rgi_ids, use_random_mb=False, tstar=tstar,
                        temp_bias=-0.5, time_scale_factors=time_scale_factors,
                        scaling_constants=scaling_const_list,
                        suffixes=suffixes, path=path, nyears=1e3)


if __name__ == '__main__':
    """ 
    """
    sensitivity()

