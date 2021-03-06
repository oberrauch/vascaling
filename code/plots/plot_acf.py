# import section
import os
import numpy as np
import pandas as pd
import xarray as xr
from statsmodels.tsa import stattools
import matplotlib.pyplot as plt
from pylab import MaxNLocator
from plots.master_colors import vas_cycle, fl_cycle

import logging

logging.basicConfig(format='%(asctime)s: %(name)s: %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    level=logging.INFO)
log = logging.getLogger('Plot ACF')


def plot_acf(ds, rgi_df, xlim=None, nlags=200, slice_start=3000,
             plot_confint=True, dir_path=None):
    """

    Parameters
    ----------
    ds
    rgi_df
    xlim
    path
    nlags
    slice_start
    plot_confint

    """
    # iterate over all above selected glaciers

    for rgi_id, glacier in rgi_df.iterrows():
        # select glacier
        rgi_id = rgi_id
        name = glacier['name']
        log.info('ACF plots for {} ({})'.format(name, rgi_id))

        # create figure and axes
        fig, ax = plt.subplots(1, 1)
        # compute acf over 1000 years
        lags = np.arange(0, nlags + 1)

        # select the complete dataset
        ds_sel = ds.sel(mb_model='random', normalized=False, rgi_id=rgi_id)
        # select time frame
        slice_end = None
        ds_sel = ds_sel.isel(time=slice(slice_start, slice_end))

        # plot zero aux line
        ax.axhline(0, c='k', ls=':')

        # iterate over all temperature biases
        for i, b in enumerate(np.sort(ds.temp_bias)):
            # get length data
            length = ds_sel.sel(temp_bias=b).length

            # FLOWLINE MODEL
            # --------------

            # compute autocorrelation and confidence intervals
            acf, confint = stattools.acf(length.sel(model='fl'), nlags=nlags,
                                         fft=True, alpha=0.01)
            # plot autocorrelation function
            ax.plot(lags[1:], acf[1:], c=fl_cycle[i],
                    label='{:+.1f} °C'.format(b))
            if plot_confint:
                # fill confidence interval
                ax.fill_between(lags[1:], confint[1:, 0] - acf[1:],
                                confint[1:, 1] - acf[1:],
                                color=fl_cycle[i], alpha=0.1)

            # V/A SCALING MODEL
            # -----------------

            # compute autocorrelation and confidence intervals
            acf, confint = stattools.acf(length.sel(model='vas'), nlags=nlags,
                                         fft=True, alpha=0.01)
            # plot autocorrelation function
            ax.plot(lags[1:], acf[1:], c=vas_cycle[i],
                    label='{:+.1f} °C'.format(b))
            if plot_confint:
                # fill confidence interval
                ax.fill_between(lags[1:], confint[1:, 0] - acf[1:],
                                confint[1:, 1] - acf[1:],
                                color=vas_cycle[i], alpha=0.1)

        # adjust axes
        if not xlim:
            xlim = [0, nlags]
        ax.set_xlim(xlim)
        ax.set_ylim([-0.55, 1.1])
        # add grid
        ax.grid()

        # get legend handles and labels
        handles, labels = ax.get_legend_handles_labels()
        title_proxy, = plt.plot(0, marker='None', linestyle='None',
                                label='dummy')

        # create list of handles and labels in correct order
        my_handles = list([title_proxy])
        my_handles.extend(handles[::2])
        my_handles.extend([title_proxy])
        my_handles.extend(handles[1::2])
        my_labels = list(["$\\bf{Flowline\ model}$"])
        my_labels.extend(labels[::2])
        my_labels.extend(["$\\bf{VAS\ model}$"])
        my_labels.extend(labels[1::2])
        # add single two-column legend
        ax.legend(my_handles, my_labels, ncol=2, loc=1)

        # labels, title, ...
        ax.set_xlabel('Lag [years]')
        ax.set_ylabel('Correlation coefficient')

        # store plot
        if not dir_path:
            dir_path = '/Users/oberrauch/work/master/plots/final_plots/acf/'
        f_name = '{}.pdf'.format(name.replace(' ', '_'))
        path = os.path.join(dir_path, f_name)
        plt.savefig(path, bbox_inches='tight')


def plot_acf_stem(ds, rgi_df, xlim=None, nlags=200, slice_start=3000,
                  plot_confint=True):
    """

    Parameters
    ----------
    ds
    rgi_df
    xlim
    path
    nlags
    slice_start
    plot_confint

    """
    # iterate over all above selected glaciers

    for rgi_id, glacier in rgi_df.iterrows():
        # select glacier
        rgi_id = rgi_id
        name = glacier['name']
        log.info('ACF plots for {} ({})'.format(name, rgi_id))

        # create figure and axes
        fig, ax = plt.subplots(1, 1)
        # compute acf over 1000 years
        lags = np.arange(0, nlags + 1)

        # select the complete dataset
        ds_sel = ds.sel(mb_model='random', normalized=False, rgi_id=rgi_id)
        # select time frame
        slice_end = None
        ds_sel = ds_sel.isel(time=slice(slice_start, slice_end))

        # plot zero aux line
        ax.axhline(0, c='k', ls=':')

        # iterate over all temperature biases
        for i, b in enumerate(np.sort(ds.temp_bias)):
            # get length data
            length = ds_sel.sel(temp_bias=b).length

            # FLOWLINE MODEL
            # --------------

            # compute autocorrelation and confidence intervals
            acf, confint = stattools.acf(length.sel(model='fl'), nlags=nlags,
                                         fft=True, alpha=0.01)
            # plot autocorrelation function
            ml, sl, bl = ax.stem(lags[1:], acf, markerfmt=f'o',
                                 linefmt=':',
                                 basefmt='None', label='{:+.1f} °C'.format(b))
            plt.setp(sl, 'color', fl_cycle[i])
            plt.setp(ml, 'color', fl_cycle[i])

            if plot_confint:
                # fill confidence interval
                ax.fill_between(lags[1:], confint[1:, 0] - acf[1:],
                                confint[1:, 1] - acf[1:],
                                color=fl_cycle[i], alpha=0.1)

            # V/A SCALING MODEL
            # -----------------

            # compute autocorrelation and confidence intervals
            acf, confint = stattools.acf(length.sel(model='vas'), nlags=nlags,
                                         fft=True, alpha=0.01)
            # plot autocorrelation function
            ml, sl, bl = ax.stem(lags[1:], acf, markerfmt=f'o',
                                 linefmt=':',
                                 basefmt='None', label='{:+.1f} °C'.format(b))
            plt.setp(sl, 'color', vas_cycle[i])
            plt.setp(ml, 'color', vas_cycle[i])

            if plot_confint:
                # fill confidence interval
                ax.fill_between(lags[1:], confint[1:, 0] - acf[1:],
                                confint[1:, 1] - acf[1:],
                                color=vas_cycle[i], alpha=0.1)

        # adjust axes
        if not xlim:
            xlim = [0, nlags]
        ax.set_xlim(xlim)
        ax.set_ylim([-1.1, 1.1])
        # add grid
        ax.grid()

        # get legend handles and labels
        handles, labels = ax.get_legend_handles_labels()
        title_proxy, = plt.plot(0, marker='None', linestyle='None',
                                label='dummy')

        # create list of handles and labels in correct order
        my_handles = list([title_proxy])
        my_handles.extend(handles[::2])
        my_handles.extend([title_proxy])
        my_handles.extend(handles[1::2])
        my_labels = list(["$\\bf{Flowline\ model}$"])
        my_labels.extend(labels[::2])
        my_labels.extend(["$\\bf{VAS\ model}$"])
        my_labels.extend(labels[1::2])
        # add single two-column legend
        ax.legend(my_handles, my_labels, ncol=2)

        # labels, title, ...
        ax.set_xlabel('Lag [years]')
        ax.set_ylabel('Correlation coefficient')

        # store plot
        dir_path = '/Users/oberrauch/work/master/plots/final_plots/acf/'
        f_name = '{}.pdf'.format(name.replace(' ', '_'))
        path = os.path.join(dir_path, f_name)
        plt.savefig(path, bbox_inches='tight')


def plot_pacf_stem(ds, rgi_df, xlim=None, path=True, nlags=200,
                   slice_start=3000,
                   plot_confint=True):
    """

    Parameters
    ----------
    ds
    rgi_df
    xlim
    path
    nlags
    slice_start
    plot_confint

    """
    # iterate over all above selected glaciers

    for rgi_id, glacier in rgi_df.iterrows():
        # select glacier
        rgi_id = rgi_id
        name = glacier['name']
        log.info('PACF plots for {} ({})'.format(name, rgi_id))

        # create figure and axes
        fig, ax = plt.subplots(1, 1)
        # compute acf over 1000 years
        lags = np.arange(0, nlags + 1)

        # select the complete dataset
        ds_sel = ds.sel(mb_model='random', normalized=False, rgi_id=rgi_id)
        # select time frame
        slice_end = None
        ds_sel = ds_sel.isel(time=slice(slice_start, slice_end))

        xoffset = 0

        # plot zero aux line
        ax.axhline(0, c='k', ls=':')

        for i, b in enumerate(np.sort(ds.temp_bias)):
            # get length data
            length = ds_sel.sel(temp_bias=b).length

            # FLOWLINE MODEL
            # --------------

            # compute autocorrelation and confidence intervals
            acf, confint = stattools.pacf(length.sel(model='fl'), nlags=nlags,
                                          alpha=0.01, method='ywmle')
            # plot autocorrelation function
            ml, sl, bl = ax.stem(lags[1:] + xoffset, acf[1:], markerfmt=f'o',
                                 linefmt=':',
                                 basefmt='None', label='{:+.1f} °C'.format(b))
            plt.setp(sl, 'color', fl_cycle[i])
            plt.setp(ml, 'color', fl_cycle[i])
            if plot_confint:
                # fill confidence interval
                ax.fill_between(lags[1:], confint[1:, 0] - acf[1:],
                                confint[1:, 1] - acf[1:],
                                color=fl_cycle[i], alpha=0.1)

            # V/A SCALING MODEL
            # -----------------

            # compute autocorrelation and confidence intervals
            acf, confint = stattools.pacf(length.sel(model='vas'), nlags=nlags,
                                          alpha=0.01, method='ywmle')
            # plot autocorrelation function
            ml, sl, bl = ax.stem(lags[1:] - xoffset, acf[1:], markerfmt=f'o',
                                 linefmt=':',
                                 basefmt='None', label='{:+.1f} °C'.format(b))
            plt.setp(sl, 'color', vas_cycle[i])
            plt.setp(ml, 'color', vas_cycle[i])
            if plot_confint:
                # fill confidence interval
                ax.fill_between(lags[1:], confint[1:, 0] - acf[1:],
                                confint[1:, 1] - acf[1:],
                                color=vas_cycle[i], alpha=0.1)

        # adjust axes
        if not xlim:
            xlim = [0, nlags]
        ax.set_xlim(xlim)
        ax.set_ylim([-1.1, 1.1])
        xa = ax.get_xaxis()
        xa.set_major_locator(MaxNLocator(integer=True))
        # add grid
        ax.grid()

        # get legend handles and labels
        handles, labels = ax.get_legend_handles_labels()
        title_proxy, = plt.plot(0, marker='None', linestyle='None',
                                label='dummy')

        # create list of handles and labels in correct order
        my_handles = list([title_proxy])
        my_handles.extend(handles[::2])
        my_handles.extend([title_proxy])
        my_handles.extend(handles[1::2])
        my_labels = list(["$\\bf{Flowline\ model}$"])
        my_labels.extend(labels[::2])
        my_labels.extend(["$\\bf{VAS\ model}$"])
        my_labels.extend(labels[1::2])
        # add single two-column legend
        ax.legend(my_handles, my_labels, ncol=2)

        # labels, title, ...
        ax.set_xlabel('Lag [years]')
        ax.set_ylabel('Correlation coefficient')

        # store plot
        dir_path = '/Users/oberrauch/work/master/plots/final_plots/pacf/'
        f_name = '{}.pdf'.format(name.replace(' ', '_'))
        path = os.path.join(dir_path, f_name)
        plt.savefig(path, bbox_inches='tight')


def plot_pacf_bars(ds, rgi_df, xlim=None, nlags=200,
                   slice_start=3000,
                   plot_confint=True):
    """

    Parameters
    ----------
    ds
    rgi_df
    xlim
    path
    nlags
    slice_start
    plot_confint

    """
    # iterate over all above selected glaciers

    for rgi_id, glacier in rgi_df.iterrows():
        # select glacier
        rgi_id = rgi_id
        name = glacier['name']
        log.info('PACF plots for {} ({})'.format(name, rgi_id))

        # create figure and axes
        fig, ax = plt.subplots(1, 1)
        # compute acf over 1000 years
        lags = np.arange(0, nlags + 1)

        # select the complete dataset
        ds_sel = ds.sel(mb_model='random', normalized=False, rgi_id=rgi_id)
        # select time frame
        slice_end = None
        ds_sel = ds_sel.isel(time=slice(slice_start, slice_end))

        # define bar width
        width = 0.4

        # plot zero aux line
        ax.axhline(0, c='k', ls=':')

        for i, b in enumerate(np.sort(ds.temp_bias)):
            # get length data
            length = ds_sel.sel(temp_bias=b).length

            # FLOWLINE MODEL
            # --------------

            # compute autocorrelation and confidence intervals
            acf, confint = stattools.pacf(length.sel(model='fl'), nlags=nlags,
                                          alpha=0.01, method='ywmle')
            # plot autocorrelation function
            ax.bar(lags - width / 2, acf, width, color=fl_cycle[i],
                   label='{:+.1f} °C'.format(b))
            if plot_confint:
                # fill confidence interval
                ax.fill_between(lags[1:], confint[1:, 0] - acf[1:],
                                confint[1:, 1] - acf[1:],
                                color=fl_cycle[i], alpha=0.1)

            # V/A SCALING MODEL
            # -----------------

            # compute autocorrelation and confidence intervals
            acf, confint = stattools.pacf(length.sel(model='vas'), nlags=nlags,
                                          alpha=0.01, method='ywmle')
            # plot autocorrelation function
            ax.bar(lags + width / 2, acf, width, color=vas_cycle[i],
                   label='{:+.1f} °C'.format(b))
            if plot_confint:
                # fill confidence interval
                ax.fill_between(lags[1:], confint[1:, 0] - acf[1:],
                                confint[1:, 1] - acf[1:],
                                color=vas_cycle[i], alpha=0.1)

        # adjust axes
        if not xlim:
            xlim = [0, nlags]
        ax.set_xlim(xlim)
        ax.set_ylim([-1.1, 1.1])
        # add grid
        ax.grid()

        # get legend handles and labels
        handles, labels = ax.get_legend_handles_labels()
        title_proxy, = plt.plot(0, marker='None', linestyle='None',
                                label='dummy')

        # create list of handles and labels in correct order
        my_handles = list([title_proxy])
        my_handles.extend(handles[::2])
        my_handles.extend([title_proxy])
        my_handles.extend(handles[1::2])
        my_labels = list(["$\\bf{Flowline\ model}$"])
        my_labels.extend(labels[::2])
        my_labels.extend(["$\\bf{VAS\ model}$"])
        my_labels.extend(labels[1::2])
        # add single two-column legend
        ax.legend(my_handles, my_labels, ncol=2)

        # labels, title, ...
        ax.set_xlabel('Lag [years]')
        ax.set_ylabel('Correlation coefficient')

        # store plot
        dir_path = '/Users/oberrauch/work/master/plots/final_plots/pacf/'
        f_name = '{}.pdf'.format(name.replace(' ', '_'))
        path = os.path.join(dir_path, f_name)
        plt.savefig(path, bbox_inches='tight')


if __name__ == '__main__':
    # specify path and read datasets
    dir_path = '/Users/oberrauch/work/master/data/' \
               + 'cluster_output/showcase_glaciers_random_climate/'
    ds = xr.open_dataset(os.path.join(dir_path, 'eq_runs.nc'))
    # sort by temperature bias
    ds = ds.sortby('temp_bias')
    # convert normalized variable from integer to boolean
    ds['normalized'] = [bool(norm) for norm in ds.normalized]

    # read showcase glaciers
    data_dir = '/Users/oberrauch/work/master/data/'
    path = os.path.join(data_dir, 'showcase_glaciers.csv')
    showcase_glaciers = pd.read_csv(os.path.join(data_dir, path), index_col=0)

    # call plotting functions
    plot_acf(ds, showcase_glaciers, nlags=500)
    plot_pacf_stem(ds, showcase_glaciers, nlags=20)
