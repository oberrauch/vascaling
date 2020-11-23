""" Compute ACF and PACF for length change signals and get order of the
autoregression and moving-average term in the data.

"""

import os
import numpy as np
import pandas as pd
import xarray as xr
from statsmodels.tsa import stattools
import matplotlib.pyplot as plt
from plots.master_colors import vas_cycle, fl_cycle

import logging

logging.basicConfig(format='%(asctime)s: %(name)s: %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    level=logging.INFO)
log = logging.getLogger('ARMA order')


def compute_arma_order_df(ds, rgi_df, nlags=200, slice_start=1000):
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
    # create an empty container
    data = list()

    # iterate over all above selected glaciers
    for rgi_id, glacier in rgi_df.iterrows():
        # create an empty container
        data_ = list()

        # select glacier
        rgi_id = rgi_id
        name = glacier['name']
        log.info('{} ({})'.format(name, rgi_id))

        # compute acf over 1000 years
        lags = np.arange(0, nlags + 1)

        # select the complete dataset
        ds_sel = ds.sel(mb_model='random', normalized=False, rgi_id=rgi_id)
        # select time frame
        slice_end = None
        ds_sel = ds_sel.isel(time=slice(slice_start, slice_end))

        # iterate over all temperature biases
        for i, b in enumerate(np.sort(ds.temp_bias)):
            # get length data
            length = ds_sel.sel(temp_bias=b).length

            # FLOWLINE MODEL
            # --------------

            # compute autocorrelation and confidence intervals
            acf, confint = stattools.acf(length.sel(model='fl'), nlags=nlags,
                                         fft=True, alpha=0.01)
            idx_significant = (acf >= confint[:, 1] - acf) | (
                        acf <= confint[:, 0] - acf)
            acf = xr.Dataset({'acf': ('lag', acf[idx_significant])},
                             coords={'lag': lags[idx_significant]})

            # compute partial autocorrelation function and confidence intervals
            pacf, confint = stattools.pacf(length.sel(model='fl'),
                                           nlags=nlags,
                                           alpha=0.01, method='ywmle')
            idx_significant = (pacf >= confint[:, 1] - pacf) | (
                        pacf <= confint[:, 0] - pacf)
            pacf = xr.Dataset({'pacf': ('lag', pacf[idx_significant])},
                              coords={'lag': lags[idx_significant]})

            ds_fl = xr.merge([acf, pacf])
            ds_fl.coords['model'] = 'fl'

            # V/A SCALING MODEL
            # -----------------

            # compute autocorrelation and confidence intervals
            acf, confint = stattools.acf(length.sel(model='vas'), nlags=nlags,
                                         fft=True, alpha=0.01)
            idx_significant = (acf >= confint[:, 1] - acf) | (
                        acf <= confint[:, 0] - acf)
            acf = xr.Dataset({'acf': ('lag', acf[idx_significant])},
                             coords={'lag': lags[idx_significant]})

            # compute partial autocorrelation function and confidence intervals
            pacf, confint = stattools.pacf(length.sel(model='vas'),
                                           nlags=nlags,
                                           alpha=0.01, method='ywmle')
            idx_significant = (pacf >= confint[:, 1] - pacf) | (
                        pacf <= confint[:, 0] - pacf)
            pacf = xr.Dataset({'pacf': ('lag', pacf[idx_significant])},
                              coords={'lag': lags[idx_significant]})

            ds_vas = xr.merge([acf, pacf])
            ds_vas.coords['model'] = 'vas'

            data__ = xr.concat([ds_fl, ds_vas], dim='model')
            data__.coords['temp_bias'] = b
            data_.append(data__)

        data_ = xr.concat(data_, dim='temp_bias')
        data_.coords['rgi_id'] = rgi_id
        data.append(data_)

    data = xr.concat(data, dim='rgi_id')

    return data


def compute_arma_order(ds, rgi_df, nlags=500, slice_start=3000):
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
    # create an empty container
    data = list()

    # iterate over all above selected glaciers
    for rgi_id, glacier in rgi_df.iterrows():
        # create an empty container
        data_ = list()

        # select glacier
        rgi_id = rgi_id
        name = glacier['name']
        log.info('{} ({})'.format(name, rgi_id))

        # compute acf over 1000 years
        lags = np.arange(0, nlags + 1)

        # select the complete dataset
        ds_sel = ds.sel(mb_model='random', normalized=False, rgi_id=rgi_id)
        # select time frame
        slice_end = None
        ds_sel = ds_sel.isel(time=slice(slice_start, slice_end))

        # iterate over all temperature biases
        for i, b in enumerate(np.sort(ds.temp_bias)):
            # get length data
            length = ds_sel.sel(temp_bias=b).length

            # FLOWLINE MODEL
            # --------------

            # compute autocorrelation and confidence intervals
            acf, confint = stattools.acf(length.sel(model='fl'), nlags=nlags,
                                         fft=True, alpha=0.01)
            idx_significant = (acf >= confint[:, 1] - acf) | (
                        acf <= confint[:, 0] - acf)
            q = lags[~idx_significant][0] - 1

            # compute partial autocorrelation function and confidence intervals
            pacf, confint = stattools.pacf(length.sel(model='fl'),
                                           nlags=nlags,
                                           alpha=0.01, method='ywmle')
            idx_significant = (pacf >= confint[:, 1] - pacf) | (
                        pacf <= confint[:, 0] - pacf)
            p = lags[~idx_significant][0] - 1
            ds_fl = xr.Dataset({'p': p, 'q': q})
            ds_fl.coords['model'] = 'fl'

            # V/A SCALING MODEL
            # -----------------

            # compute autocorrelation and confidence intervals
            acf, confint = stattools.acf(length.sel(model='vas'), nlags=nlags,
                                         fft=True, alpha=0.01)
            idx_significant = (acf >= confint[:, 1] - acf) | (
                    acf <= confint[:, 0] - acf)
            q = lags[~idx_significant][0]

            # compute partial autocorrelation function and confidence intervals
            pacf, confint = stattools.pacf(length.sel(model='vas'),
                                           nlags=nlags,
                                           alpha=0.01, method='ywmle')
            idx_significant = (pacf >= confint[:, 1] - pacf) | (
                    pacf <= confint[:, 0] - pacf)
            p = lags[~idx_significant][0]
            ds_vas = xr.Dataset({'p': p, 'q': q})
            ds_vas.coords['model'] = 'vas'

            data__ = xr.concat([ds_fl, ds_vas], dim='model')
            data__.coords['temp_bias'] = b
            data_.append(data__)

        data_ = xr.concat(data_, dim='temp_bias')
        data_.coords['rgi_id'] = rgi_id
        data.append(data_)

    data = xr.concat(data, dim='rgi_id')

    return data


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
    f_name = 'showcase_glaciers.csv'
    showcase_glaciers = pd.read_csv(os.path.join(data_dir, f_name),
                                    index_col=0)

    data = compute_arma_order(ds, showcase_glaciers)
    f_name = 'arma_order.nc'
    data.to_netcdf(os.path.join(data_dir, f_name))
