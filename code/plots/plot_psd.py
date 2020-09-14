# import section
import numpy as np
import pandas as pd
import xarray as xr
from scipy import signal
import matplotlib.pyplot as plt

import logging
logging.basicConfig(format='%(asctime)s: %(name)s: %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    level=logging.INFO)
log = logging.getLogger('plot psd')

# specify path and read datasets
path = '/Users/oberrauch/work/master/data/random_climate_same_tstar/eq_runs.nc'
ds = xr.open_dataset(path)
# sort by temperature bias
ds = ds.sortby('temp_bias')
ds['normalized'] = [bool(norm) for norm in ds.normalized]

# mass balance data set
path = '/Users/oberrauch/work/master/data/random_climate_same_tstar/mb_output.nc'
ds_mb = xr.open_dataset(path)
# sort by temperature bias
ds_mb = ds_mb.sortby('temp_bias')

# read showcase glaciers
showcase_glaciers = pd.read_csv('/Users/oberrauch/work/master/data/showcase_glaciers.csv', index_col=0)

# define color cycles
vas_cycle = np.array(['#f7ca18', '#f39c12', '#c0392b', '#22313f', '#4d13d1', '#59abe3'])
fl_cycle = np.array(["#4ecdc4", "#1f3a93", "#a537fd", "#26a65b", "#00e640", "#0093ac"])


def plot_psd():
    # iterate over all above selected glaciers
    for rgi_id, glacier in showcase_glaciers.iterrows():
        # select glacier
        rgi_id = rgi_id
        name = glacier['name']
        # create figure and axes
        fig, ax = plt.subplots(1, 1, figsize=[8, 6])
        # select from complete dataset
        ds_sel = ds.sel(normalized=False,
                        mb_model='random',
                        rgi_id=rgi_id)
        # truncate spinup if necessary
        ds_sel = ds_sel.isel(time=slice(0, None))

        # plots for the FLOWLINE model
        for i, b in enumerate(np.sort(ds.temp_bias)):
            # select values by temperature bias
            ds_tmp = ds_sel.sel(temp_bias=b).length
            # compute the power of the signel per frequency band
            sig = ds_tmp.sel(model='fl').values.flatten()
            freqs, psd = signal.welch(sig)
            # convert frequency into periods
            periods = 1 / freqs
            ax.loglog(periods, psd, label='{:+.1f} °C'.format(b), c=fl_cycle[i], lw=2)

        # plots for the V/A SCALING model
        for i, b in enumerate(np.sort(ds.temp_bias)):
            # select values by temperature bias
            ds_tmp = ds_sel.sel(temp_bias=b).length
            # compute the power of the signel per frequency band
            sig = ds_tmp.sel(model='vas').values.flatten()
            freqs, psd = signal.welch(sig)
            # convert frequency into periods
            periods = 1 / freqs
            ax.loglog(periods, psd, label='{:+.1f} °C'.format(b), c=vas_cycle[i], lw=2)

        # get legend handles and labels
        handles, labels = ax.get_legend_handles_labels()
        title_proxy, = plt.plot(0, marker='None', linestyle='None', label='dummy')

        # add a seperate legend for each model
        # leg_fl = ax.legend(handles[:3], labels[:3], bbox_to_anchor=(1, 1), loc='upper left')
        # leg_fl.set_title('Flowline model', {'weight': 'bold'})
        # leg_vas = ax.legend(handles[3:], labels[3:], bbox_to_anchor=(1, 0), loc='lower left')
        # leg_vas.set_title('VAS model', {'weight': 'bold'})
        # ax.add_artist(leg_fl)

        # create list of handles and labels in correct order
        my_handles = list([title_proxy])
        my_handles.extend(handles[:3])
        my_handles.extend([title_proxy])
        my_handles.extend(handles[3:])
        my_labels = list(["$\\bf{Flowline\ model}$"])
        my_labels.extend(labels[:3])
        my_labels.extend(["$\\bf{VAS\ model}$"])
        my_labels.extend(labels[3:])
        # add single two-column legend
        ax.legend(my_handles, my_labels, ncol=2)

        # add grid
        ax.grid(which='both')
        # invert x-axis
        ax.invert_xaxis()

        # title, labels, ...
        ax.set_title('PSD of {} length changes under random climate'.format(name),
                     {'weight': 'bold'})
        ax.set_xlabel('Period [years]')
        ax.set_ylabel('Power density [m$^2$/yr]')

        # store plot
        f_path = '/Users/oberrauch/work/master/plots/final_plots/psd/{}.pdf'.format(name.replace(' ', '_'))
        plt.savefig(f_path, bbox_inches='tight')


def plot_psd_length_mb():
    # iterate over all above selected glaciers
    for rgi_id, glacier in showcase_glaciers.iterrows():
        # select glacier
        rgi_id = rgi_id
        name = glacier['name']
        # create figure and axes
        fig, ax = plt.subplots(1, 2, figsize=[10, 6])
        # select from complete dataset
        ds_sel = ds.sel(normalized=False,
                        mb_model='random',
                        rgi_id=rgi_id)
        ds_mb_sel = ds_mb.sel(rgi_id=rgi_id)
        # truncate spinup if necessary
        ds_sel = ds_sel.isel(time=slice(0, None))
        ds_mb_sel = ds_mb_sel.isel(year=slice(0, None))

        # plot spectral density of length changes in left subplot
        # -------------------------------------------------------

        # plots for the FLOWLINE model
        for i, b in enumerate(np.sort(ds.temp_bias)):
            # select values by temperature bias
            ds_tmp = ds_sel.sel(temp_bias=b).length

            # compute the power of the signel per frequency band
            sig = ds_tmp.sel(model='fl').values.flatten()
            freqs, psd = signal.welch(sig)
            ax[0].loglog(1 / freqs, psd, label='{:+.1f} °C'.format(b), c=fl_cycle[i], lw=2)
            ax[1].loglog(1 / freqs, psd, label='{:+.1f} °C'.format(b), c='lightgrey', lw=1)

        # plots for the V/A SCALING model
        for i, b in enumerate(np.sort(ds.temp_bias)):
            # select values by temperature bias
            ds_tmp = ds_sel.sel(temp_bias=b).length

            # compute the power of the signel per frequency band
            sig = ds_tmp.sel(model='vas').values.flatten()
            freqs, psd = signal.welch(sig)
            ax[0].loglog(1 / freqs, psd, label='{:+.1f} °C'.format(b), c=vas_cycle[i], lw=2)
            ax[1].loglog(1 / freqs, psd, label='{:+.1f} °C'.format(b), c='lightgrey', lw=1)

        # plot spectral density of mass balance in right subplot
        # ------------------------------------------------------
        for j, b in enumerate(np.sort(ds_mb.temp_bias)):
            # select values by temperature bias
            ds_tmp = ds_mb_sel.sel(temp_bias=b).spec_mb

            # compute the power of the signel per frequency band
            sig = ds_tmp.sel(model='fl').values.flatten()
            freqs, psd = signal.welch(sig)
            ax[1].loglog(1 / freqs, psd, label='{:+.1f} °C'.format(b), c=fl_cycle[j], lw=2)
            # compute the power of the signel per frequency band
            sig = ds_tmp.sel(model='vas').values.flatten()
            freqs, psd = signal.welch(sig)
            ax[1].loglog(1 / freqs, psd, label='{:+.1f} °C'.format(b), c=vas_cycle[j], lw=2)


                # get legend handles and labels
        handles, labels = ax[0].get_legend_handles_labels()
        title_proxy, = plt.plot(0, marker='None', linestyle='None', label='dummy')

        # add a seperate legend for each model
        # leg_fl = ax.legend(handles[:3], labels[:3], bbox_to_anchor=(1, 1), loc='upper left')
        # leg_fl.set_title('Flowline model', {'weight': 'bold'})
        # leg_vas = ax.legend(handles[3:], labels[3:], bbox_to_anchor=(1, 0), loc='lower left')
        # leg_vas.set_title('VAS model', {'weight': 'bold'})
        # ax.add_artist(leg_fl)

        # create list of handles and labels in correct order
        my_handles = list([title_proxy])
        my_handles.extend(handles[:3])
        my_handles.extend([title_proxy])
        my_handles.extend(handles[3:])
        my_labels = list(["$\\bf{Flowline\ model}$"])
        my_labels.extend(labels[:3])
        my_labels.extend(["$\\bf{VAS\ model}$"])
        my_labels.extend(labels[3:])
        # add single two-column legend
        ax[1].legend(my_handles, my_labels, ncol=2, loc='lower center')

        # add grid
        ax[0].grid(which='both')
        ax[1].grid(which='both')

        # invert x-axis
        ax[0].invert_xaxis()
        ax[1].invert_xaxis()

        # title, labels, ...
        ax[0].set_title('Glacier length', {'weight': 'bold'})
        ax[1].set_title('Specific mass balance', {'weight': 'bold'})
        ax[0].set_xlabel('Frequency [year$^{-1}$]')
        ax[1].set_xlabel('Frequency [year$^{-1}$]')
        ax[0].set_ylabel('Power density [m$^2$/yr]')
        ax[0].set_ylabel('Power density [(mm w.e./yr)$^2$/yr]')

        # set figure title
        fig.suptitle('PSD - {} under a random climate scenario'.format(name), weight='bold', fontsize=18)

        # store plot
        f_path = '/Users/oberrauch/work/master/plots/final_plots/psd/{}_mb.pdf'.format(name.replace(' ', '_'))
        plt.savefig(f_path, bbox_inches='tight')


if __name__ == '__main__':
    plot_psd()
    plot_psd_length_mb()
