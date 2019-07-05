# Equilibrium runs
# --------------------
# This script runs the VAS model for a single (or more) glacier(s) with a
# constant or random massbalance model, performing equilibrium experiments.
# An additional comparison with the OGGM flowline model is performed.
# TODO

# import externals libraries
import os
import shutil
import matplotlib.pyplot as plt

# import the needed OGGM modules
from oggm import cfg, utils, workflow
from oggm.core import gis, climate, flowline, vascaling


def normalize_ds_with_start(ds):
    """ Normalize all data variables with their respective first entry.
    Return a new xarray.Dataset.

    Parameters
    ----------
    ds: xarray.Dataset

    Returns
    -------

    """
    # copy the given dataset
    ds_ = ds.copy()
    # iterate over all data variables
    for var in ds_.data_vars:
        # add normalize variable as attribute
        ds_.attrs[var + '_0'] = ds_[var][0]
        # normalize with start value
        ds_[var] = ds_[var] / ds_[var][0]
    # return normalized dataset
    return ds_


def equilibrium_run_vas(rgi_ids, use_random_mb=True, use_mean=True,
                        plot=False, figure_title='', paths=True,  **kwargs):
    """

    Returns
    -------

    """
    # specify suffixes
    suffixes = ['_normal', '_bias_p', '_bias_n', '_longtime']

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
    shutil.rmtree(wdir)
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
    workflow.execute_entity_task(vascaling.local_t_star, gdirs)

    if use_random_mb:
        # RANDOM MASS BALANCE MODEL
        # -------------------------

        # use t* as center year, even if specified differently
        kwargs['y0'] = None
        # run for 500 years if not specified otherwise
        kwargs.setdefault('nyears', 500)
        # set random seed to get reproducible results
        kwargs.setdefault('seed', 12)

        # run RandomMassBalance model centered around t*, once without
        # temperature bias and once with positive and negative temperature bias
        # of 0.5 °C each.
        workflow.execute_entity_task(vascaling.run_random_climate, gdirs,
                                     **kwargs, output_filesuffix=suffixes[0])
        workflow.execute_entity_task(vascaling.run_random_climate, gdirs,
                                     **kwargs, temperature_bias=+0.5,
                                     output_filesuffix=suffixes[1])
        workflow.execute_entity_task(vascaling.run_random_climate, gdirs,
                                     **kwargs, temperature_bias=-0.5,
                                     output_filesuffix=suffixes[2])
        # run over a longer time period
        kwargs['nyears'] = 3000
        workflow.execute_entity_task(vascaling.run_random_climate, gdirs,
                                     **kwargs, output_filesuffix=suffixes[3])
    else:
        # CONSTANT MASS BALANCE MODEL
        # -------------------------

        # use t* as center year, even if specified differently
        kwargs['y0'] = None
        # run for 500 years if not specified otherwise
        kwargs.setdefault('nyears', 500)

        # the bias is defined to be zero during the calibration process,
        # which is why we don't use it here to reproduce the results
        cfg.PARAMS['use_bias_for_run'] = False
        # run RandomMassBalance model centered around t*, once without
        # temperature bias and once with positive and negative temperature bias
        # of 0.5 °C each.
        workflow.execute_entity_task(vascaling.run_constant_climate, gdirs,
                                     **kwargs, output_filesuffix=suffixes[0])
        workflow.execute_entity_task(vascaling.run_constant_climate, gdirs,
                                     **kwargs, temperature_bias=+0.5,
                                     output_filesuffix=suffixes[1])
        workflow.execute_entity_task(vascaling.run_constant_climate, gdirs,
                                     **kwargs, temperature_bias=-0.5,
                                     output_filesuffix=suffixes[2])
        # run over a longer time period
        kwargs['nyears'] = 3000
        workflow.execute_entity_task(vascaling.run_constant_climate, gdirs,
                                     **kwargs, output_filesuffix=suffixes[3])

    # compile the run output and normalize glacier geometries with start value
    ds_normal = utils.compile_run_output(gdirs, filesuffix=suffixes[0])
    ds_bias_p = utils.compile_run_output(gdirs, filesuffix=suffixes[1])
    ds_bias_n = utils.compile_run_output(gdirs, filesuffix=suffixes[2])
    ds_longtime = utils.compile_run_output(gdirs, filesuffix=suffixes[3])

    # normalize
    if use_mean:
        ds_normal = normalize_ds_with_start(ds_normal).mean(dim='rgi_id')
        ds_bias_p = normalize_ds_with_start(ds_bias_p).mean(dim='rgi_id')
        ds_bias_n = normalize_ds_with_start(ds_bias_n).mean(dim='rgi_id')
        ds_longtime = normalize_ds_with_start(ds_longtime).mean(dim='rgi_id')
    else:
        ds_normal = normalize_ds_with_start(ds_normal.sum(dim='rgi_id'))
        ds_bias_p = normalize_ds_with_start(ds_bias_p.sum(dim='rgi_id'))
        ds_bias_n = normalize_ds_with_start(ds_bias_n.sum(dim='rgi_id'))
        ds_longtime = normalize_ds_with_start(ds_longtime.sum(dim='rgi_id'))

    if plot:
        # create figure and axes
        fig, [ax0, ax1, ax2] = plt.subplots(3, 1, figsize=[6, 8])
        # plot the evolution of glacier volume
        ax0.plot(ds_normal.volume, label='equilibrium', c='#2e3131')
        ax0.plot(ds_bias_p.volume, label='+0.5 °C', c='#f22613')
        ax0.plot(ds_bias_n.volume, label='-0.5 °C', c='#1f3a93')
        ax0.axhline(ds_normal.volume[0], c='k', ls=':', lw=0.8, label='initial value')
        ax0.set_xticklabels('')
        ax0.set_ylabel('Relative volume')
        ax0.set_title(figure_title)
        ax0.legend(bbox_to_anchor=(0.5, 0), loc=9, ncol=4)

        # plot the evolution of glacier area
        ax1.plot(ds_normal.area, label='equilibrium', c='#2e3131')
        ax1.plot(ds_bias_p.area, label='+0.5 °C', c='#f22613')
        ax1.plot(ds_bias_n.area, label='-0.5 °C', c='#1f3a93')
        ax1.axhline(ds_normal.area[0], c='k', ls=':', lw=0.8, label='initial value')
        ax1.set_xticklabels('')
        ax1.set_ylabel('Relative area')
        ax1.legend(bbox_to_anchor=(0.5, 0), loc=9, ncol=4)

        # plot the evolution of glacier length
        ax2.plot(ds_normal.length, label='equilibrium', c='#2e3131')
        ax2.plot(ds_bias_p.length, label='+0.5 °C', c='#f22613')
        ax2.plot(ds_bias_n.length, label='-0.5 °C', c='#1f3a93')
        ax2.axhline(ds_normal.length[0], c='k', ls=':', lw=0.8, label='initial value')
        ax2.set_xlabel('Years of evolution')
        ax2.set_ylabel('Relative length')

        # show plot and store to file
        plt.show()
        fig.savefig('vas_relative_geometries.pdf', bbox_inches='tight')

    # store normalized datasets
    if paths:
        if paths is True:
            paths = [os.path.join(cfg.PATHS['working_dir'], 'normalized_'
                                  + suffix + '.nc')
                     for suffix in suffixes]
        for path, suffix in zip(paths, suffixes):
            ds = eval('ds' + suffix)
            ds.to_netcdf(path)

    return {'normal': ds_normal, 'bias_p': ds_bias_p, 'bias_n': ds_bias_n,
            'longtime': ds_longtime}


def equilibrium_run_fl(rgi_ids, use_random_mb=True, use_mean=True,
                       plot=False, figure_title='', paths=True, **kwargs):
    """

    Returns
    -------

    """
    # specify suffixes
    suffixes = ['_normal', '_bias_p', '_bias_n', '_longtime']

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
    shutil.rmtree(wdir)
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

    # run gis tasks
    workflow.gis_prepro_tasks(gdirs)
    # run climate tasks
    workflow.climate_tasks(gdirs)
    # run inversion tasks
    workflow.inversion_tasks(gdirs)
    # finalize preprocessing
    workflow.execute_entity_task(flowline.init_present_time_glacier, gdirs)

    if use_random_mb:
        # RANDOM MASS BALANCE MODEL
        # -------------------------

        # use t* as center year, even if specified differently
        kwargs['y0'] = None
        # run for 500 years if not specified otherwise
        kwargs.setdefault('nyears', 500)
        # set random seed to get reproducible results
        kwargs.setdefault('seed', 12)

        # run RandomMassBalance model centered around t*, once without
        # temperature bias and once with positive and negative temperature bias
        # of 0.5 °C each.
        workflow.execute_entity_task(flowline.run_random_climate, gdirs,
                                     **kwargs, output_filesuffix=suffixes[0])
        workflow.execute_entity_task(flowline.run_random_climate, gdirs,
                                     **kwargs, temperature_bias=+0.5,
                                     output_filesuffix=suffixes[1])
        workflow.execute_entity_task(flowline.run_random_climate, gdirs,
                                     **kwargs, temperature_bias=-0.5,
                                     output_filesuffix=suffixes[2])
        # run over a longer time period
        kwargs['nyears'] = 3000
        workflow.execute_entity_task(flowline.run_random_climate, gdirs,
                                     **kwargs, output_filesuffix=suffixes[3])
    else:
        # CONSTANT MASS BALANCE MODEL
        # -------------------------

        # use t* as center year, even if specified differently
        kwargs['y0'] = None
        # run for 500 years if not specified otherwise
        kwargs.setdefault('nyears', 500)

        # the bias is defined to be zero during the calibration process,
        # which is why we don't use it here to reproduce the results
        cfg.PARAMS['use_bias_for_run'] = False
        # run RandomMassBalance model centered around t*, once without
        # temperature bias and once with positive and negative temperature bias
        # of 0.5 °C each.
        workflow.execute_entity_task(flowline.run_constant_climate, gdirs,
                                     **kwargs, output_filesuffix=suffixes[0])
        workflow.execute_entity_task(flowline.run_constant_climate, gdirs,
                                     **kwargs, temperature_bias=+0.5,
                                     output_filesuffix=suffixes[1])
        workflow.execute_entity_task(flowline.run_constant_climate, gdirs,
                                     **kwargs, temperature_bias=-0.5,
                                     output_filesuffix=suffixes[2])
        # run over a longer time period
        kwargs['nyears'] = 3000
        workflow.execute_entity_task(flowline.run_constant_climate, gdirs,
                                     **kwargs, output_filesuffix=suffixes[3])

    # compile the run output and normalize glacier geometries with start value
    ds_normal = utils.compile_run_output(gdirs, filesuffix=suffixes[0])
    ds_bias_p = utils.compile_run_output(gdirs, filesuffix=suffixes[1])
    ds_bias_n = utils.compile_run_output(gdirs, filesuffix=suffixes[2])
    ds_longtime = utils.compile_run_output(gdirs, filesuffix=suffixes[3])

    # normalize
    if use_mean:
        ds_normal = normalize_ds_with_start(ds_normal).mean(dim='rgi_id')
        ds_bias_p = normalize_ds_with_start(ds_bias_p).mean(dim='rgi_id')
        ds_bias_n = normalize_ds_with_start(ds_bias_n).mean(dim='rgi_id')
        ds_longtime = normalize_ds_with_start(ds_longtime).mean(dim='rgi_id')
    else:
        ds_normal = normalize_ds_with_start(ds_normal.sum(dim='rgi_id'))
        ds_bias_p = normalize_ds_with_start(ds_bias_p.sum(dim='rgi_id'))
        ds_bias_n = normalize_ds_with_start(ds_bias_n.sum(dim='rgi_id'))
        ds_longtime = normalize_ds_with_start(ds_longtime.sum(dim='rgi_id'))

    if plot:
        # create figure and axes
        fig, [ax0, ax1, ax2] = plt.subplots(3, 1, figsize=[6, 8])
        # plot the evolution of glacier volume
        ax0.plot(ds_normal.volume, label='equilibrium', c='#2e3131')
        ax0.plot(ds_bias_p.volume, label='+0.5 °C', c='#f22613')
        ax0.plot(ds_bias_n.volume, label='-0.5 °C', c='#1f3a93')
        ax0.axhline(ds_normal.volume[0], c='k', ls=':', lw=0.8, label='initial value')
        ax0.set_xticklabels('')
        ax0.set_ylabel('Relative volume')
        ax0.set_title(figure_title)
        ax0.legend(bbox_to_anchor=(0.5, 0), loc=9, ncol=4)

        # plot the evolution of glacier area
        ax1.plot(ds_normal.area, label='equilibrium', c='#2e3131')
        ax1.plot(ds_bias_p.area, label='+0.5 °C', c='#f22613')
        ax1.plot(ds_bias_n.area, label='-0.5 °C', c='#1f3a93')
        ax1.axhline(ds_normal.area[0], c='k', ls=':', lw=0.8, label='initial value')
        ax1.set_xticklabels('')
        ax1.set_ylabel('Relative area')
        ax1.legend(bbox_to_anchor=(0.5, 0), loc=9, ncol=4)

        # plot the evolution of glacier length
        ax2.plot(ds_normal.length, label='equilibrium', c='#2e3131')
        ax2.plot(ds_bias_p.length, label='+0.5 °C', c='#f22613')
        ax2.plot(ds_bias_n.length, label='-0.5 °C', c='#1f3a93')
        ax2.axhline(ds_normal.length[0], c='k', ls=':', lw=0.8, label='initial value')
        ax2.set_xlabel('Years of evolution')
        ax2.set_ylabel('Relative length')

        # show plot and store to file
        plt.show()
        fig.savefig('fl_relative_geometries.pdf', bbox_inches='tight')

    # store normalized datasets
    if paths:
        if paths is True:
            paths = [os.path.join(cfg.PATHS['working_dir'], 'normalized'
                                  + suffix + '.nc')
                     for suffix in suffixes]
        for path, suffix in zip(paths, suffixes):
            ds = eval('ds' + suffix)
            ds.to_netcdf(path)

    return {'normal': ds_normal, 'bias_p': ds_bias_p, 'bias_n': ds_bias_n,
            'longtime': ds_longtime}


if __name__ == '__main__':
    # define RGI IDs
    rgi_ids = ['RGI60-11.00897', 'RGI60-11.01270']
    # perform equilibrium experiments
    vas_ds = equilibrium_run_vas(rgi_ids)
    fl_ds = equilibrium_run_fl(rgi_ids)
