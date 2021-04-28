# import section
import os
import numpy as np
import pandas as pd
import xarray as xr
from scipy import signal
import matplotlib.pyplot as plt

from plots.master_colors import vas_cycle, fl_cycle

import logging
logging.basicConfig(format='%(asctime)s: %(name)s: %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    level=logging.INFO)
log = logging.getLogger('Plot PSD')


def plot_psd(ds, rgi_df, plot_periods=False, spinup_time=0, norm=False):
    # iterate over all above selected glaciers
    for rgi_id, glacier in rgi_df.iterrows():
        # select glacier
        rgi_id = rgi_id
        name = glacier['name']
        log.info(f'Plotting PSD for {name}')
        # create figure and axes
        fig, ax = plt.subplots(1, 1, figsize=[8, 6])
        # select from complete dataset
        ds_sel = ds.sel(normalized=False,
                        mb_model='random',
                        rgi_id=rgi_id)
        # truncate spinup if necessary
        ds_sel = ds_sel.isel(time=slice(int(spinup_time), None))

        # plots for the FLOWLINE model
        for i, b in enumerate(np.sort(ds.temp_bias)):
            # select values by temperature bias
            ds_tmp = ds_sel.sel(temp_bias=b).length
            # compute the power of the signal per frequency band
            sig = ds_tmp.sel(model='fl').values.flatten()
            nperseg = int(len(sig) / 10)
            freqs, psd = signal.welch(sig, nperseg=nperseg)
            # convert frequency into periods
            if plot_periods:
                freqs = 1 / freqs
            # normalize power density
            if norm:
                psd = psd / psd[0]
            ax.loglog(freqs, psd, label='{:+.1f} °C'.format(b), c=fl_cycle[i],
                      lw=2)

        # plots for the V/A SCALING model
        for i, b in enumerate(np.sort(ds.temp_bias)):
            # select values by temperature bias
            ds_tmp = ds_sel.sel(temp_bias=b).length
            # compute the power of the signel per frequency band
            sig = ds_tmp.sel(model='vas').values.flatten()
            freqs, psd = signal.welch(sig, nperseg=nperseg)
            # convert frequency into periods
            if plot_periods:
                freqs = 1 / freqs
            # normalize power density
            if norm:
                psd = psd / psd[0]
            ax.loglog(freqs, psd, label='{:+.1f} °C'.format(b), c=vas_cycle[i],
                      lw=2)

        # make x-axis tight
        ax.autoscale(enable=True, axis='x', tight=True)

        # get legend handles and labels
        handles, labels = ax.get_legend_handles_labels()
        title_proxy, = plt.plot(0, marker='None', linestyle='None',
                                label='dummy')

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
        ax.legend(my_handles, my_labels, ncol=2, loc=3)

        # add grid
        ax.grid(which='both')
        # invert x-axis
        if plot_periods:
            ax.invert_xaxis()

        # title, labels, ...
        if plot_periods:
            ax.set_xlabel('Period [yr]')
        else:
            ax.set_xlabel('Frequency [yr$^{-1}$]')
        if norm:
            ax.set_ylabel('Normalized power density')
        else:
            ax.set_ylabel('Power density [m$^2$/yr$^{-1}$]')

        # store plot
        dir_path = '/Users/oberrauch/work/master/plots/final_plots/psd/'
        f_name = '{}{}.pdf'.format(name.replace(' ', '_'),
                                   '_norm' if norm else '')
        # plt.savefig(os.path.join(dir_path, f_name), bbox_inches='tight')
        print(ax.get_ylim())


def plot_psd_mb(ds_mb, rgi_df, plot_periods=False, spinup_time=0, norm=False):
    # iterate over all above selected glaciers
    for rgi_id, glacier in rgi_df.iterrows():
        # select glacier
        rgi_id = rgi_id
        name = glacier['name']
        log.info(f'Plotting PSD for {name}')
        # create figure and axes
        fig, ax = plt.subplots(1, 1, figsize=[8, 6])
        # select from complete dataset
        ds_sel = ds_mb.sel(rgi_id=rgi_id)
        # truncate spinup if necessary
        ds_sel = ds_sel.isel(year=slice(int(spinup_time), None))

        # plots for the FLOWLINE model
        for i, b in enumerate(np.sort(ds_mb.temp_bias)):
            # select values by temperature bias
            ds_tmp = ds_sel.sel(temp_bias=b).spec_mb
            # compute the power of the signal per frequency band
            sig = ds_tmp.sel(model='fl').values.flatten()
            nperseg = int(len(sig) / 10)
            freqs, psd = signal.welch(sig, nperseg=nperseg)
            # convert frequency into periods
            if plot_periods:
                freqs = 1 / freqs
            # normalize power density
            if norm:
                psd = psd / psd[0]
            ax.loglog(freqs, psd, label='{:+.1f} °C'.format(b), c=fl_cycle[i],
                      lw=2)

        # plots for the V/A SCALING model
        for i, b in enumerate(np.sort(ds_mb.temp_bias)):
            # select values by temperature bias
            ds_tmp = ds_sel.sel(temp_bias=b).spec_mb
            # compute the power of the signel per frequency band
            sig = ds_tmp.sel(model='vas').values.flatten()
            freqs, psd = signal.welch(sig, nperseg=nperseg)
            # convert frequency into periods
            if plot_periods:
                freqs = 1 / freqs
            # normalize power density
            if norm:
                psd = psd / psd[0]
            ax.loglog(freqs, psd, label='{:+.1f} °C'.format(b), c=vas_cycle[i],
                      lw=2)

        # make x-axis tight
        ax.autoscale(enable=True, axis='x', tight=True)

        # get legend handles and labels
        handles, labels = ax.get_legend_handles_labels()
        title_proxy, = plt.plot(0, marker='None', linestyle='None',
                                label='dummy')

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
        ax.legend(my_handles, my_labels, ncol=2, loc=3)

        # add grid
        ax.grid(which='both')

        ax.set_ylim((0.010493453462322543, 282327512.32164896))
        # invert x-axis
        if plot_periods:
            ax.invert_xaxis()

        # title, labels, ...
        if plot_periods:
            ax.set_xlabel('Period [yr]')
        else:
            ax.set_xlabel('Frequency [yr$^{-1}$]')
        if norm:
            ax.set_ylabel('Normalized power density')
        else:
            ax.set_ylabel('Power density [(mm w.e./yr)$^2$/yr$^{-1}$]')

        # store plot
        dir_path = '/Users/oberrauch/work/master/plots/final_plots/psd/'
        f_name = '{}{}_mb.pdf'.format(name.replace(' ', '_'),
                                   '_norm' if norm else '')
        plt.savefig(os.path.join(dir_path, f_name), bbox_inches='tight')


def plot_psd_length_mb(ds, ds_mb, rgi_df, plot_periods=False, spinup_time=0):
    # iterate over all above selected glaciers
    for rgi_id, glacier in rgi_df.iterrows():
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
        ds_sel = ds_sel.isel(time=slice(int(spinup_time), None))
        ds_mb_sel = ds_mb_sel.isel(year=slice(int(spinup_time), None))

        # plot spectral density of length changes in left subplot
        # -------------------------------------------------------

        # plots for the FLOWLINE model
        for i, b in enumerate(np.sort(ds.temp_bias)):
            # select values by temperature bias
            ds_tmp = ds_sel.sel(temp_bias=b).length

            # compute the power of the signel per frequency band
            sig = ds_tmp.sel(model='fl').values.flatten()
            nperseg = int(len(sig) / 10)
            freqs, psd = signal.welch(sig, nperseg=nperseg)
            # convert frequencies into period
            if plot_periods:
                freqs = 1 / freqs
            ax[0].loglog(freqs, psd, label='{:+.1f} °C'.format(b),
                         c=fl_cycle[i], lw=2)
            ax[1].loglog(freqs, psd, label='{:+.1f} °C'.format(b),
                         c='lightgrey', lw=1)

        # plots for the V/A SCALING model
        for i, b in enumerate(np.sort(ds.temp_bias)):
            # select values by temperature bias
            ds_tmp = ds_sel.sel(temp_bias=b).length

            # compute the power of the signel per frequency band
            sig = ds_tmp.sel(model='vas').values.flatten()
            freqs, psd = signal.welch(sig, nperseg=nperseg)
            # convert frequencies into period
            if plot_periods:
                freqs = 1 / freqs
            ax[0].loglog(freqs, psd, label='{:+.1f} °C'.format(b),
                         c=vas_cycle[i], lw=2)
            ax[1].loglog(freqs, psd, label='{:+.1f} °C'.format(b),
                         c='lightgrey', lw=1)

        # plot spectral density of mass balance in right subplot
        # ------------------------------------------------------
        for j, b in enumerate(np.sort(ds_mb.temp_bias)):
            # select values by temperature bias
            ds_tmp = ds_mb_sel.sel(temp_bias=b).spec_mb

            # compute the power of the signel per frequency band
            sig = ds_tmp.sel(model='fl').values.flatten()
            nperseg = int(len(sig) / 10)
            freqs, psd = signal.welch(sig, nperseg=nperseg)
            # convert frequencies into period
            if plot_periods:
                freqs = 1 / freqs
            ax[1].loglog(freqs, psd, label='{:+.1f} °C'.format(b),
                         c=fl_cycle[j], lw=2)
            # compute the power of the signel per frequency band
            sig = ds_tmp.sel(model='vas').values.flatten()
            freqs, psd = signal.welch(sig, nperseg=nperseg)
            # convert frequencies into period
            if plot_periods:
                freqs = 1 / freqs
            ax[1].loglog(freqs, psd, label='{:+.1f} °C'.format(b),
                         c=vas_cycle[j], lw=2)

        # get legend handles and labels
        handles, labels = ax[0].get_legend_handles_labels()
        title_proxy, = plt.plot(0, marker='None', linestyle='None',
                                label='dummy')

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

        # make x-axis tight
        [ax_.autoscale(enable=True, axis='x', tight=True) for ax_ in ax]

        # invert x-axis
        if plot_periods:
            ax[0].invert_xaxis()
            ax[1].invert_xaxis()

        # title, labels, ...
        ax[0].set_title('Glacier length', {'weight': 'bold'})
        ax[1].set_title('Specific mass balance', {'weight': 'bold'})
        if plot_periods:
            ax[0].set_xlabel('Period [yr]')
            ax[1].set_xlabel('Period [yr]')
        else:
            ax[0].set_xlabel('Frequency [yr$^{-1}$]')
            ax[1].set_xlabel('Frequency [yr$^{-1}$]')
        ax[0].set_ylabel('Power density [m$^2$/yr$^{-1}$]')
        ax[0].set_ylabel('Power density [(mm w.e./yr)$^2$/yr$^{-1}$]')

        # set figure title
        fig.suptitle('PSD - {} under a random climate scenario'.format(name),
                     weight='bold', fontsize=18)

        # store plot
        f_path = '/Users/oberrauch/work/master/plots/final_plots/psd/{}_mb.pdf'.format(
            name.replace(' ', '_'))
        plt.savefig(f_path, bbox_inches='tight')


if __name__ == '__main__':
    # Configure logger to display time and script name in a nice way
    logging.basicConfig(format='%(asctime)s: %(name)s: %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        level=logging.INFO)
    log = logging.getLogger('plot psd')

    # specify path and read datasets
    data_dir = '/Users/oberrauch/work/master/data/' \
               + 'cluster_output/showcase_glaciers_random_climate/'
    path = os.path.join(data_dir, 'eq_runs.nc')
    ds = xr.load_dataset(path)
    data_dir = '/Users/oberrauch/work/master/data/' \
               + 'cluster_output/showcase_glaciers_random_climate/'
    path = os.path.join(data_dir, 'mb_output.nc')
    ds_mb = xr.load_dataset(path)

    # convert normalized variable from integer to boolean
    ds['normalized'] = [bool(norm) for norm in ds.normalized]

    # sort by temperature bias
    ds = ds.sortby('temp_bias')
    ds_mb = ds_mb.sortby('temp_bias')

    # read showcase glaciers
    data_dir = '/Users/oberrauch/work/master/data/'
    path = os.path.join(data_dir, 'showcase_glaciers.csv')
    showcase_glaciers = pd.read_csv(os.path.join(data_dir, path), index_col=0)

    # define spinup time
    spinup_time = 3e3

    # call plotting functions
    # plot_psd(ds, showcase_glaciers, spinup_time=spinup_time, norm=False)
    # plot_psd_length_mb(ds, ds_mb, showcase_glaciers, spinup_time)
    plot_psd_mb(ds_mb, showcase_glaciers, spinup_time=spinup_time)
