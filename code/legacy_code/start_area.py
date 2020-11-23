""" Find model start areas (for 1851), using the iterative process from Marzeion, et.
al. (2012). The routine includes all needed steps, only the RGI ID and the
glacier's name must be specified.

The script runs the task for all 'reference' glaciers.
"""

## Import section
# import externals libs
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# import the needed OGGM modules
import oggm
from oggm import cfg, utils, workflow
from oggm.core import gis, climate
import oggm_vas as vascaling


def seek_start_area(rgi_id, name, show=False, path='', ref=np.NaN,
                    adjust_term_elev=False, legend=True,
                    instant_geometry_change=False):
    """ Set up an VAS model from scratch and run/test the start area seeking
    tasks. The result is a plot showing the modeled glacier area evolution for
    different start values. The plots can be displayed and/or stored to file.

    Parameters
    ----------
    rgi_id: string
        RGI ID denoting the glacier on which to perform the tasks
    name: string
        Name og glacier, since it is not always given (or correct) in RGI
    show: bool, optional, default=False
        Flag deciding whether or not to show the created plots.
    path: string, optional, default=''
        Path under which the modeled area plot should be stored.
    ref: float, optional, default=np.NaN
        Historic (1851) reference area with which a reference model run is
        performed.

    """
    # Initialization and load default parameter file
    vascaling.initialize()

    # compute RGI region and version from RGI IDs
    # assuming they all are all the same
    rgi_region = (rgi_id.split('-')[-1]).split('.')[0]
    rgi_version = (rgi_id.split('-')[0])[-2:-1]

    # specify working directory and output directory
    working_dir = os.path.abspath('../working_directories/start_area/')
    # output_dir = os.path.abspath('./vas_run_output')
    output_dir = os.path.abspath('../data/vas_run_output')
    # create working directory
    utils.mkdir(working_dir, reset=False)
    utils.mkdir(output_dir)
    # set path to working directory
    cfg.PATHS['working_dir'] = working_dir
    # set RGI version and region
    cfg.PARAMS['rgi_version'] = rgi_version
    # define how many grid points to use around the glacier,
    # if you expect the glacier to grow large use a larger border
    cfg.PARAMS['border'] = 20
    # we use HistAlp climate data
    cfg.PARAMS['baseline_climate'] = 'HISTALP'
    # set the mb hyper parameters accordingly
    cfg.PARAMS['prcp_scaling_factor'] = 1.75
    cfg.PARAMS['temp_melt'] = -1.75
    cfg.PARAMS['run_mb_calibration'] = False

    # the bias is defined to be zero during the calibration process,
    # which is why we don't use it here to reproduce the results
    cfg.PARAMS['use_bias_for_run'] = True

    # get/downlaod the rgi entity including the outline shapefile
    rgi_df = utils.get_rgi_glacier_entities([rgi_id])
    # set name, if not delivered with RGI
    if rgi_df.loc[int(rgi_id[-5:]) - 1, 'Name'] is None:
        rgi_df.loc[int(rgi_id[-5:]) - 1, 'Name'] = name

    # get and set path to intersect shapefile
    intersects_db = utils.get_rgi_intersects_region_file(region=rgi_region)
    cfg.set_intersects_db(intersects_db)

    # initialize the GlacierDirectory
    gdir = workflow.init_glacier_directories(rgi_df)[0]

    # # DEM and GIS tasks
    # # get the path to the DEM file (will download if necessary)
    # dem = utils.get_topo_file(gdir.cenlon, gdir.cenlat)
    # print('DEM source: {}, path to DEM file: {}'.format(dem[1], dem[0][0]))
    # # set path in config file
    # cfg.PATHS['dem_file'] = dem[0][0]
    # cfg.PARAMS['border'] = 10
    # cfg.PARAMS['use_intersects'] = False

    # run GIS tasks
    gis.define_glacier_region(gdir)
    gis.glacier_masks(gdir)

    # process climate data
    climate.process_climate_data(gdir)

    #  compute local t* and the corresponding mu*
    vascaling.local_t_star(gdir)

    # create mass balance model
    mb_mod = vascaling.VAScalingMassBalance(gdir)

    # look at specific mass balance over climate data period
    min_hgt, max_hgt = vascaling.get_min_max_elevation(gdir)
    y0 = 1851
    y1 = 2014

    # run scalar minimization
    minimize_res = vascaling.find_start_area(gdir,
                                             adjust_term_elev=adjust_term_elev,
                                             instant_geometry_change=instant_geometry_change)
    # print(minimize_res)

    # stop script if minimization was not successful
    if minimize_res.status and False:
        sys.exit(minimize_res.status)

    # instance glacier with today's values
    model_ref = vascaling.VAScalingModel(year_0=gdir.rgi_date,
                                         area_m2_0=gdir.rgi_area_m2,
                                         min_hgt=min_hgt, max_hgt=max_hgt,
                                         mb_model=mb_mod)

    # instance guessed starting areas
    num = 9
    area_guess = np.linspace(1e6, np.floor(gdir.rgi_area_m2 * 2), num,
                             endpoint=True)
    # create empty containers
    area_list = list()
    volume_list = list()
    spec_mb_list = list()

    # iterate over all starting areas
    for area_ in area_guess:
        # instance iteration model
        model_guess = vascaling.VAScalingModel(year_0=gdir.rgi_date,
                                               area_m2_0=gdir.rgi_area_m2,
                                               min_hgt=min_hgt,
                                               max_hgt=max_hgt,
                                               mb_model=mb_mod)
        # set new starting values
        model_guess.create_start_glacier(area_, y0,
                                         adjust_term_elev=adjust_term_elev)
        # run model and save years and area
        best_guess_ds = model_guess.run_until_and_store(
            year_end=model_ref.year,
            instant_geometry_change=instant_geometry_change)
        # create series and store in container
        area_list.append(best_guess_ds.area_m2.to_dataframe()['area_m2'])
        volume_list.append(best_guess_ds.volume_m3.to_dataframe()['volume_m3'])
        spec_mb_list.append(best_guess_ds.spec_mb.to_dataframe()['spec_mb'])

    # create DataFrame
    area_df = pd.DataFrame(area_list, index=['{:.2f}'.format(a / 1e6)
                                             for a in area_guess])
    area_df.index.name = 'Start Area [km$^2$]'

    volume_df = pd.DataFrame(volume_list, index=['{:.2f}'.format(a / 1e6)
                                                 for a in area_guess])
    volume_df.index.name = 'Start Area [km$^2$]'

    # set up model with resulted starting area
    model = vascaling.VAScalingModel(year_0=model_ref.year_0,
                                     area_m2_0=model_ref.area_m2_0,
                                     min_hgt=model_ref.min_hgt_0,
                                     max_hgt=model_ref.max_hgt,
                                     mb_model=model_ref.mb_model)
    model.create_start_glacier(minimize_res.x, y0,
                               adjust_term_elev=adjust_term_elev)

    # run model with best guess initial area
    best_guess_ds = model.run_until_and_store(year_end=model_ref.year,
                                              instant_geometry_change=instant_geometry_change)
    # run model with historic reference area
    if ref:
        model.reset()
        model.create_start_glacier(ref * 1e6, y0,
                                   adjust_term_elev=adjust_term_elev)
        ref_ds = model.run_until_and_store(year_end=model_ref.year,
                                           instant_geometry_change=instant_geometry_change)

    # create figure and add axes
    fig = plt.figure(figsize=[5, 5])
    ax = fig.add_axes([0.125, 0.075, 0.85, 0.9])

    # plot model output
    ax = (area_df / 1e6).T.plot(legend=False, colormap='Spectral', ax=ax)

    # plot best guess
    ax.plot(best_guess_ds.time, best_guess_ds.area_m2 / 1e6, color='k',
            ls='--', lw=1.2,
            label=f'{best_guess_ds.area_m2.isel(time=0).values/1e6:.2f} km$^2$ (best result)')
    # plot reference
    if ref:
        ax.plot(ref_ds.time, ref_ds.area_m2 / 1e6, color='k',
                ls='-.', lw=1.2,
                label=f'{ref_ds.area_m2.isel(time=0).values/1e6:.2f} km$^2$ (1850 ref.)')

    # plot 2003 reference line
    ax.axhline(model_ref.area_m2_0 / 1e6, c='k',
               ls=':',
               label=f'{model_ref.area_m2_0/1e6:.2f} km$^2$ ({gdir.rgi_date} obs.)')

    # add legend
    if legend:
        handels, labels = ax.get_legend_handles_labels()
        labels[:-3] = [r'{} km$^2$'.format(l) for l in labels[:-3]]
        leg = ax.legend(handels, labels, loc='upper right', ncol=2)
        # leg.set_title('Start area $A_0$', prop={'size': 12})

    # replot best guess estimate and reference (in case it lies below another
    # guess)
    ax.plot(best_guess_ds.time, best_guess_ds.area_m2 / 1e6, color='k',
            ls='--', lw=1.2)
    if ref:
        ax.plot(ref_ds.time, ref_ds.area_m2 / 1e6, color='k', ls='-.', lw=1.2)

    # labels, title
    ax.set_xlim([best_guess_ds.time.values[0], best_guess_ds.time.values[-1]])
    ax.set_xlabel('')
    ax.set_ylabel('Glacier area [km$^2$]')

    # save figure to file
    if path:
        fig.savefig(path)

    # show plot
    if show:
        plt.show()
    plt.clf()

    # plot and store volume evolution
    (volume_df / 1e9).T.plot(legend=False, colormap='viridis')
    plt.gcf().savefig(path[:-4] + '_volume.pdf')


if __name__ == '__main__':
    """ The script runs the task for all 'reference' glaciers. """
    # define list of glaciers
    rgi_ids = ['RGI60-11.00897', 'RGI60-11.00787', 'RGI60-11.01450']
    names = ['Hintereisferner', 'Kesselwandferner', 'Großer Aletschgletscher']
    # historic area from RGI/GLIMS
    ref_areas_1850 = [15.4, 5.11, 105.61]

    for rgi_id, name, ref_area in zip(rgi_ids, names, ref_areas_1850):
        # define path for plots
        path = '/Users/oberrauch/work/master/plots/start_area/{}_{}.pdf'

        # define file name
        f_path = path.format(rgi_id, 'original')
        # original implementation
        seek_start_area(rgi_id, name, path=f_path,
                        ref=ref_area, show=False,
                        adjust_term_elev=False, legend=False)

        # define file name
        f_path = path.format(rgi_id, 'overturn')
        # physical consistent implementation
        seek_start_area(rgi_id, name, path=f_path,
                        ref=ref_area, show=False,
                        adjust_term_elev=True, legend=True)

        # define file name
        f_path = path.format(rgi_id, 'instant')
        # physical consistent implementation with instant geometry update
        seek_start_area(rgi_id, name, path=f_path,
                        ref=ref_area, show=False, instant_geometry_change=True,
                        adjust_term_elev=True, legend=False)
