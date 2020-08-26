# import section
import numpy as np
import pandas as pd
import xarray as xr
from statsmodels.tsa import stattools
import matplotlib.pyplot as plt

import logging
logging.basicConfig(format='%(asctime)s: %(name)s: %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    level=logging.INFO)
log = logging.getLogger('plot acf')

# specify path and read datasets
path = '/Users/oberrauch/work/master/data/random_climate_same_tstar/eq_runs.nc'
ds = xr.open_dataset(path)
# sort by temperature bias
ds = ds.sortby('temp_bias')
ds['normalized'] = [bool(norm) for norm in ds.normalized]

# define color cycles
vas_cycle = np.array(['#f7ca18', '#f39c12', '#c0392b', '#22313f', '#4d13d1', '#59abe3'])
fl_cycle = np.array(["#4ecdc4", "#1f3a93", "#a537fd", "#26a65b", "#00e640", "#0093ac"])

# iterate over all above selected glaciers
showcase_glaciers = pd.read_csv('/Users/oberrauch/work/master/data/showcase_glaciers.csv', index_col=0)
for rgi_id, glacier in showcase_glaciers.iterrows():
    # select glacier
    rgi_id = rgi_id
    name = glacier['name']
    log.info('ACF plots for {} ({})'.format(name, rgi_id))

    # create figure and axes
    fig, ax = plt.subplots(1, 1)
    # compute acf over 1000 years
    nlags = 1000

    # select the complete dataset
    ds_sel = ds.sel(mb_model='random',
                    normalized=False,
                    rgi_id=rgi_id)

    for i, b in enumerate(np.sort(ds.temp_bias)):
        # get length data
        length = ds_sel.sel(temp_bias=b).length
        # plot autocorrelation
        ax.plot(stattools.acf(length.sel(model='fl'), nlags=nlags, fft=True),
                c=fl_cycle[i], label='{:+.1f} °C'.format(b))
        ax.plot(stattools.acf(length.sel(model='vas'), nlags=nlags, fft=True),
                c=vas_cycle[i], label='{:+.1f} °C'.format(b))

    # aux line
    ax.axhline(0, c='k', ls=':')
    # adjust axes
    ax.set_xlim([0, nlags])
    ylim = ax.get_ylim()
    ax.set_ylim([min(ylim), 1])
    # add grid
    ax.grid()

    # get legend handles and labels
    handles, labels = ax.get_legend_handles_labels()
    # add a seperate legend for each model
    leg_fl = ax.legend(handles[::2], labels[::2], bbox_to_anchor=(1, 1), loc='upper left')
    leg_fl.set_title('Flowline model', {'weight': 'bold'})
    leg_vas = ax.legend(handles[1::2], labels[1::2], bbox_to_anchor=(1, 0), loc='lower left')
    leg_vas.set_title('VAS model', {'weight': 'bold'})
    ax.add_artist(leg_fl)

    # labels, title, ...
    ax.set_title('ACF of {} length under random climate'.format(name), {'weight': 'bold'})
    ax.set_xlabel('Lag [years]')
    ax.set_ylabel('Correlation coefficient')

    # store plot
    f_path = '/Users/oberrauch/work/master/plots/final_plots/acf/{}.pdf'.format(name.replace(' ', '_'))
    plt.savefig(f_path, bbox_inches='tight')
