# import section
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

import logging
logging.basicConfig(format='%(asctime)s: %(name)s: %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    level=logging.INFO)
log = logging.getLogger('plot timeseries')


def plot_time_series(ds, var, title='', suptitle='', xlim=None):
    # define color cycles
    # define color cycles
    vas_cycle = np.array(['#f7ca18', '#f39c12', '#c0392b', '#22313f', '#4d13d1', '#59abe3'])
    fl_cycle = np.array(["#4ecdc4", "#1f3a93", "#a537fd", "#26a65b", "#00e640", "#0093ac"])

    # plot relative volume change
    fig, [ax0, ax1] = plt.subplots(1, 2, figsize=[10, 5])

    # flowline model
    ds.sel(model='fl')[var].plot(hue='temp_bias', ax=ax0, add_legend=False, color='lightgray', lw=0.5)
    # vas model
    ax0.set_prop_cycle('color', vas_cycle)
    handles_vas = ds.sel(model='vas')[var].plot(hue='temp_bias', ax=ax0, add_legend=False, lw=2)
    labels_vas = ['{:+.1f} °C'.format(bias) for bias in ds.temp_bias.values]

    # set axes limits
    if xlim:
        ax0.set_xlim(xlim)
    else:
        ax0.set_xlim([ds.time.min(), ds.time.max()])
    ylim = ax0.get_ylim()
    ax0.set_ylim([ylim[0] * 0.5, ylim[1]])

    # title, labels, legend
    ax0.set_title('Volume/area scaling model')
    ax0.set_xlabel('Years of model evolution')
    ax0.legend(handles_vas, labels_vas, title='Temperature bias',
               bbox_to_anchor=(0.5, 0), loc=8, ncol=3)
    if ds.normalized:
        # add ylabel
        ax0.set_ylabel('Relative {}'.format(var))
        # aux line
        ax0.axhline(1, lw=0.8, ls=':', c='k')
    else:
        # add ylabel
        unit = 'm' if var == 'length' else 'm$^2$' if var == 'area' else 'm$^3$'
        ax0.set_ylabel('Glacier {} [{}]'.format(var, unit))
        # aux line
        ax0.axhline(ds.sel(model='vas')[var].isel(time=0).mean(), lw=0.8, ls=':', c='k')

    ax0.grid()

    # vas model
    ds.sel(model='vas')[var].plot(hue='temp_bias', ax=ax1, add_legend=False, color='lightgray', lw=0.5)
    # flowline model
    ax1.set_prop_cycle('color', fl_cycle)
    handles_fl = ds.sel(model='fl')[var].plot(hue='temp_bias', ax=ax1, add_legend=False)
    labels_fl = ['{:+.1f} °C'.format(bias) for bias in ds.temp_bias.values]

    # set axes limits
    if xlim:
        ax1.set_xlim(xlim)
    else:
        ax1.set_xlim([ds.time.min(), ds.time.max()])
    ylim = ax1.get_ylim()
    ax1.set_ylim([ylim[0] * 0.5, ylim[1]])

    # title, labels, legend
    ax1.set_title('Flowline model')
    ax1.set_xlabel('Years of model evolution')
    ax1.legend(handles_fl, labels_fl, title='Temperature bias',
               bbox_to_anchor=(0.5, 0), loc=8, ncol=3)
    ax1.yaxis.tick_right()
    ax1.yaxis.set_label_position('right')
    if ds.normalized:
        # add ylabel
        ax1.set_ylabel('Relative {}'.format(var))
        # aux line
        ax1.axhline(1, lw=0.8, ls=':', c='k')
    else:
        # add ylabel
        unit = 'm' if var == 'length' else 'm$^2$' if var == 'area' else 'm$^3$'
        ax1.set_ylabel('Glacier {} [{}]'.format(var, unit))
        # aux line
        ax1.axhline(ds.sel(model='fl')[var].isel(time=0).mean(), lw=0.8, ls=':', c='k')

    ax1.grid()

    # add suptitle
    fig.suptitle(suptitle, fontsize=15)


def plot_random_single_glaciers():
    # specify path and read datasets
    path = '/Users/oberrauch/work/master/data/random_climate_same_tstar/eq_runs.nc'
    ds = xr.open_dataset(path)
    # sort by temperature bias
    ds = ds.sortby('temp_bias')
    ds['normalized'] = [bool(norm) for norm in ds.normalized]

    # define glaciers of interest
    ref_glaciers = pd.read_csv('/Users/oberrauch/work/master/data/showcase_glaciers.csv', index_col=0)

    # iterate over all above selected glaciers
    for rgi_id, glacier in ref_glaciers.iterrows():
        # select glacier
        rgi_id = rgi_id
        name = glacier['name']
        log.info('Plots for {} ({})'.format(name, rgi_id))

        # iterate over all selected variales
        variables = ['length', 'volume']
        for var in variables:
            for norm in [True, False]:
                log.info('Plotting {} {} time series'
                         .format('normalized' if norm else 'absolute',
                                 var))
                suptitle = '{} evolution under random climate'.format(name)
                plot_time_series(ds.sel(mb_model='random',
                                        normalized=norm,
                                        rgi_id=rgi_id),
                                 var=var, suptitle=suptitle, xlim=[0, 1e3])
                f_path = '/Users/oberrauch/work/master/plots/final_plots/time_series/single_glaciers/'
                f_path += '{}_{}_{}.pdf'.format(var, 'norm' if norm else 'abs',
                                                name.replace(' ', '_'))
                plt.savefig(f_path, bbox_inches='tight')
                plt.close(plt.gcf())

