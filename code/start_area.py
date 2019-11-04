# Run and find start ares
# Set up an OGGM/VAS run from scratch and test the start area seeking tasks.
# The only thing to specify ist the RGI ID and the glacier's name, the rest
# should run without any adjustments...

## Import section
# import externals libs
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# import the needed OGGM modules
import oggm
from oggm import cfg, utils
from oggm.utils import get_rgi_glacier_entities
from oggm.core import gis, climate
from oggm.core import vascaling


def seek_start_area(rgi_id, name, show=False, path='', ref=np.NaN):
    """

    :param rgi_id:
    :param name:
    :param show:
    :param path:
    :param ref:
    :return:
    """
    ## Initialization
    # load parameter file
    cfg.initialize()

    # get/downlaod the rgi entity including the outline shapefile
    rgi_df = get_rgi_glacier_entities([rgi_id])
    # set name, if not delivered with RGI
    # if rgi_df.loc[int(rgi_id[-5:])-1, 'Name'] is None:
    rgi_df.loc[int(rgi_id[-5:])-1, 'Name'] = name

    # select single entry
    rgi_entity = rgi_df.iloc[0]
    # visualize
    rgi_df.plot()
    plt.title(rgi_entity.Name)
    # show plot
    if show:
        plt.show()
    plt.clf()


    ## Glacier Directory
    # specify the working directory and define the glacier directory
    cfg.PATHS['working_dir'] = './'
    gdir = oggm.GlacierDirectory(rgi_entity, reset=True)

    ## DEM and GIS tasks
    # get the path to the DEM file (will download if necessary)
    dem = utils.get_topo_file(gdir.cenlon, gdir.cenlat)
    print('DEM source: {}, path to DEM file: {}'.format(dem[1], dem[0][0]))

    # set path in config file
    cfg.PATHS['dem_file'] = dem[0][0]
    cfg.PARAMS['border'] = 10
    cfg.PARAMS['use_intersects'] = False

    # run GIS tasks
    gis.define_glacier_region(gdir, entity=rgi_entity)
    gis.glacier_masks(gdir)

    ## Climate data
    # set mb calibration parameters
    cfg.PARAMS['baseline_climate'] = 'HISTALP'
    cfg.PARAMS['prcp_scaling_factor'] = 1.75
    cfg.PARAMS['temp_melt'] = -1.75
    # process HistAlp climate data
    climate.process_histalp_data(gdir)

    ## Mass balance model
    #  compute local t* and the corresponding mu*
    vascaling.local_t_star(gdir)
    # see calibration results
    print(gdir.read_json('vascaling_mustar'))

    # create mass balance model
    mb_mod = vascaling.VAScalingMassBalance(gdir)

    # look at specific mass balance over climate data period
    min_hgt, max_hgt = vascaling.get_min_max_elevation(gdir)
    y0 = 1802
    y1 = 2014
    years = np.arange(y0, y1)
    mb = list()
    for y in years:
        mb.append(mb_mod.get_specific_mb(min_hgt, max_hgt, y))

    # visualize
    plt.plot(years, mb)
    plt.axhline(0, c='k', ls=':', lw=0.8)
    plt.title('Modeled mass balance - {}'.format(name))
    plt.ylabel('Mass balance [mm w.e. yr$^{-1}$]')
    # show plot
    if show:
        plt.show()
    plt.clf()

    ## Find start area
    # run scalar minimization
    minimize_res = vascaling.find_start_area(gdir)
    print(minimize_res)

    # stop script if minimization was not successful
    if minimize_res.status:
        sys.exit(minimize_res.status)

    # instance glacier with today's values
    model_ref = vascaling.VAScalingModel(year_0=gdir.rgi_date,
                                         area_m2_0=gdir.rgi_area_m2,
                                         min_hgt=min_hgt, max_hgt=max_hgt,
                                         mb_model=mb_mod)

    # instance guessed starting areas
    num = 15
    area_guess = np.linspace(100, gdir.rgi_area_m2*2,  num, endpoint=True)
    # create empty containers
    iteration_list = list()
    spec_mb_list = list()

    # iterate over all starting areas
    for area_ in area_guess:
        # instance iteration model
        model_guess = vascaling.VAScalingModel(year_0=gdir.rgi_date,
                                               area_m2_0=gdir.rgi_area_m2,
                                               min_hgt=min_hgt, max_hgt=max_hgt,
                                               mb_model=mb_mod)
        # set new starting values
        model_guess.create_start_glacier(area_, 1851)
        # run model and save years and area
        best_guess_ds = model_guess.run_until_and_store(year_end=model_ref.year)
        # create series and store in container
        iteration_list.append(best_guess_ds.area_m2.to_dataframe()['area_m2'])
        spec_mb_list.append(best_guess_ds.spec_mb.to_dataframe()['spec_mb'])

    # create DataFrame
    iteration_df = pd.DataFrame(iteration_list, index=['{:.2f}'.format(a/1e6) for a in area_guess])
    iteration_df.index.name = 'Start Area [km$^2$]'

    # set up model with resulted starting area
    model = vascaling.VAScalingModel(year_0=model_ref.year_0,
                                     area_m2_0=model_ref.area_m2_0,
                                     min_hgt=model_ref.min_hgt_0,
                                     max_hgt=model_ref.max_hgt,
                                     mb_model=model_ref.mb_model)
    model.create_start_glacier(minimize_res.x, 1851)

    # run model with best guess initial area
    best_guess_ds = model.run_until_and_store(year_end=model_ref.year)
    # run model with historic reference area
    if ref:
        model.reset()
        model.create_start_glacier(ref*1e6, 1851)
        ref_ds = model.run_until_and_store(year_end=model_ref.year)

    # create figure and add axes
    fig = plt.figure(figsize=[11, 6])
    ax = fig.add_axes([0.1, 0.1, 0.65, 0.8])
    # plot best guess
    ax.plot(best_guess_ds.time, best_guess_ds.area_m2 / 1e6, color='k',
            ls='--', lw=1.2, label='$A_0$ with best result')
    # plot reference
    if ref:
        ax.plot(ref_ds.time, ref_ds.area_m2 / 1e6, color='k',
                ls='-.', lw=1.2, label='ref. area from 1850')
    # plot model output
    ax = (iteration_df / 1e6).T.plot(legend=False, colormap='Spectral', ax=ax)
    # plot 2003 reference line
    ax.axhline(model_ref.area_m2_0 / 1e6, c='k',
               ls=':', label='measured area in {}'.format(gdir.rgi_date))
    # add legend
    handels, labels = ax.get_legend_handles_labels()
    labels[2:-1] = [r'{} km$^2$'.format(l) for l in labels[2:-1]]
    leg = ax.legend(handels, labels, bbox_to_anchor=(1.025, 0.5), loc='center left')
    leg.set_title('Start area $A_0$', prop={'size': 12})
    leg._legend_box.align = 'left'

    # replot best guess estimate and reference (in case it lies below another guess)
    ax.plot(best_guess_ds.time, best_guess_ds.area_m2 / 1e6, color='k',
            ls='--', lw=1.2)
    if ref:
        ax.plot(ref_ds.time, ref_ds.area_m2 / 1e6, color='k', ls='-.', lw=1.2)

    # labels, title
    ax.set_xlim([best_guess_ds.time.values[0], best_guess_ds.time.values[-1]])
    ax.set_xlabel('')
    ax.set_ylabel('Glacier area [km$^2$]')
    fig.suptitle('Modeled glacier area - {}'.format(rgi_entity.Name), fontsize='14')

    # save figure to file
    if path:
        fig.savefig(path)

    # show plot
    if show:
        plt.show()
    plt.clf()


if __name__ == '__main__':
    # get list of demo glaciers and select those in the Alps/HistAlp domain
    demo_glaciers = pd.read_csv('/Users/oberrauch/oggm-fork/oggm/data/demo_glaciers.csv', index_col=1)
    demo_glaciers = demo_glaciers[demo_glaciers.RGIId.str.contains('11.')]
    # add historic area from RGI/GLIMS
    demo_glaciers['area_1850_km2'] = [15.4, 5.11, 105.61, np.NaN, np.NaN, 10.12, 33.41]
    demo_glaciers
    for _, glacier in demo_glaciers.iterrows():
        fn = glacier.Name.lower().replace(' ', '_')
        fn = '/Users/oberrauch/work/master/plots/start_area/{}.pdf'.format(fn)
        seek_start_area(glacier.RGIId, glacier.Name, path=fn,
                        ref=glacier.area_1850_km2)
        break

