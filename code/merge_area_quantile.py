import pandas as pd
import os, glob
import shutil
import numpy as np
from scipy import stats
import xarray as xr
import matplotlib.pyplot as plt
import progressbar


def mkdir(path, reset=False):
    """Checks if directory exists and if not, create one.

    Parameters
    ----------
    reset: erase the content of the directory if exists

    Returns
    -------
    the path
    """

    if reset and os.path.exists(path):
        shutil.rmtree(path)
        # deleting stuff takes time
        while os.path.exists(path):  # check if it still exists
            pass
    try:
        os.makedirs(path)
    except FileExistsError:
        pass
    return path


def extended_vas(f1, f2, sel_ids=None):
    """A function which merges together the historical runs with the forward
    runs for convenience.

    Parameters:
    -----------
    f1: str, path to historical run *.nc file
    f2: str, path to forward run *.nc file
    sel_ids: string array-like, optional, default=None
        RGI IDs to exclude from datasets

    Returns:
    --------
    the combined xr.DataSet

    """

    # open both *.nc files
    with xr.open_dataset(f1) as ds1, xr.open_dataset(f2) as ds2:
        # load and subset both datasets
        ds1 = ds1.load()
        ds1 = ds1.isel(rgi_id=~ds1.rgi_id.isin(sel_ids))
        ds2 = ds2.load()
        ds2 = ds2.isel(rgi_id=~ds2.rgi_id.isin(sel_ids))

        # some housekeeping, getting rid of suffixes and certain variables
        ds1 = ds1.rename_vars(
            dict((k, k.replace('_ext', '')) for k in ds1.data_vars))

        # subset historical run
        ds1 = ds1.sel(time=slice(2000, 2019))

        # combine both datasets into one along the time dimension
        ds = xr.concat([ds1, ds2], dim='time')

        # ajust rate variables
        dt_vars = [v for v in ds.data_vars if '_dt' in v]
        for vn in dt_vars:
            ds[vn].data[1:, :] = ds[vn].data[:-1, :]
            ds[vn].data[0, :] = np.NaN

        # subset for 21st century
        ds = ds.sel(time=slice(2000, None))

    return ds


def extended_fl(f1, f2, sel_ids=None):
    """A function which merges together the historical runs with the
    forward runs for convenience.

    Parameters:
    -----------
    f1: str, path to historical run *.nc file
    f2: str, path to forward run *.nc file
    sel_ids: string array-like, optional, default=None
        RGI IDs to exclude from datasets

    Returns:
    --------
    the combined xr.DataSet

    """

    # open both *.nc files
    with xr.open_dataset(f1) as ds1, xr.open_dataset(f2) as ds2:
        # load and subset both datasets
        ds1 = ds1.load()
        ds1 = ds1.isel(rgi_id=~ds1.rgi_id.isin(sel_ids))
        ds2 = ds2.load()
        ds2 = ds2.isel(rgi_id=~ds2.rgi_id.isin(sel_ids))

        # some housekeeping, getting rid of suffixes and certain variables
        ds1 = ds1.rename_vars(
            dict((k, k.replace('_ext', '')) for k in ds1.data_vars))
        ds1 = ds1.drop_vars('volume_fixed_geom')

        # compute calving flux and add as variable
        vn = 'calving_dt'
        for ds in [ds1, ds2]:
            ds[vn] = ds['calving'].copy(deep=True)
            ds[vn].attrs['description'] = 'Yearly calving flux (past year)'
            ds[vn].attrs['unit'] = 'm3 yr-1'
            ds[vn].data[:-1, :] = ds[vn].data[1:, :] - ds[vn].data[0:-1, :]
            ds[vn].data[-1, :] = np.NaN

        # compute fluxes for volume and volume below sea level
        # add as variables to dataset
        for vno in ['volume', 'volume_bsl']:
            vn = vno + '_dt'
            for ds in [ds1, ds2]:
                ds[vn] = ds[vno].copy(deep=True)
                ds[vn].attrs['description'] += 'change (past year)'
                ds[vn].attrs['unit'] += ' yr-1'
                ds[vn].data[:-1, :] = ds[vn].data[1:, :] - ds[vn].data[0:-1, :]
                ds[vn].data[-1, :] = np.NaN

        # adjsut calving with last value
        ds1['calving'] = ds1['calving'] - ds1['calving'].isel(time=-1)

        # subset historical run
        ds1 = ds1.sel(time=slice(2000, 2019))

        # combine both datasets into one along the time dimension
        ds = xr.concat([ds1, ds2], dim='time')

        # ajust rate variables
        dt_vars = [v for v in ds.data_vars if '_dt' in v]
        for vn in dt_vars:
            ds[vn].data[1:, :] = ds[vn].data[:-1, :]
            ds[vn].data[0, :] = np.NaN

        # subset for 21st century
        ds = ds.sel(time=slice(2000, None))

    return ds


def merge_vas():
    # specify region
    rgi_reg = int(os.environ['RGI_REG'])
    ref_year = 2020

    # define input directory
    idir = '/home/users/moberrauch/cmip6_run/cmip6_output/'

    # get list of failed glaciers
    failed = pd.read_csv('/home/users/moberrauch/cmip6_run/'
                         'failed_glaciers_vas.csv',
                         index_col=0)
    failed = failed.values.reshape(len(failed))

    # specify path to output *.nc files for historical runs
    f_past_tpl = '/home/www/moberrauch/prepro/run_output_{}.nc'

    # specify dummy path to glacier statistics
    stas_dirpath = '/home/www/oggm/gdirs/oggm_v1.4/exps/thesis_vas/' \
                   'RGI62/b_080/L5/summary/glacier_statistics_{}.csv'

    glacier_stats = pd.read_csv(stas_dirpath.format(f'{rgi_reg:02d}'),
                                index_col=0, low_memory=False)

    # some logging info
    print(f'RGI region: {rgi_reg}', flush=True)

    # relative volume of given region, in realtion to reference year
    out_frac_small = dict()
    out_frac_large = dict()
    # absolute volume of given region
    out_abs_small = dict()
    out_abs_large = dict()

    # get RGI IDs for large glaciers
    quantile = 0.95
    area_threshold = glacier_stats.rgi_area_km2.quantile(quantile)
    small_ids = glacier_stats[
        (glacier_stats.rgi_area_km2 <= area_threshold)].index.values
    large_ids = glacier_stats[
        (glacier_stats.rgi_area_km2 > area_threshold)].index.values

    # make sure the RGI ID is a string
    rgi_reg = '{:02d}'.format(rgi_reg)
    # get RGI IDs for failed glaciers (see section above)
    sel_ids = [r for r in failed if rgi_reg + '.' in r]

    sel_large_ids = np.unique(np.concatenate([sel_ids, small_ids]))
    sel_small_ids = np.unique(np.concatenate([sel_ids, large_ids]))

    # get path to historic run for given RGI region
    f_past = f_past_tpl.format(rgi_reg)

    # define paths and create output directories
    odir_frac = '/home/users/moberrauch/run_output/' \
                'merge_cru_cmip_vas/frac_output/' + rgi_reg + '/'
    odir_abs = '/home/users/moberrauch/run_output/' \
               'merge_cru_cmip_vas/abs_output/' + rgi_reg + '/'
    mkdir(odir_frac)
    mkdir(odir_abs)

    # get path to RGI region and the *.nc files
    # for all the GCMs and SSPs
    rdir = idir + 'RGI' + rgi_reg + '/'
    all_ncs = glob.glob(rdir + '*/*.nc')

    # select certain ssps
    ncs = list()
    for ssp_ in ['ssp126', 'ssp245', 'ssp370', 'ssp585']:
        ncs.extend([nc for nc in all_ncs if ssp_ in nc])

    # iterate over all *.nc files, i.e. all combinations of region/GCM/SSP
    for nc in progressbar.progressbar(sorted(ncs)):

        # get SSP and GCM description
        ssp = os.path.basename(nc).split('_')[1].replace('.nc', '')
        gcm = os.path.basename(nc).split('_')[0]
        # define key for dataset
        key_both = ssp + '_' + gcm
        # create new dataframe for give SSP, if it doesn't already exist
        # and add it to nested dictionary for absolute and relative volume
        if ssp not in out_frac_small:
            out_frac_small[ssp] = pd.DataFrame()
            out_frac_large[ssp] = pd.DataFrame()
            out_abs_small[ssp] = pd.DataFrame()
            out_abs_large[ssp] = pd.DataFrame()

        # merge historical run and forward run
        ds_small = extended_vas(f_past, nc, sel_ids=sel_small_ids)
        ds_large = extended_vas(f_past, nc, sel_ids=sel_large_ids)

        # sum volume over all glaciers and convert into pd.Series
        all_small_vols = ds_small.volume.sum(dim='rgi_id').to_series()
        all_large_vols = ds_large.volume.sum(dim='rgi_id').to_series()

        # use longest time index for regional volume dictionaries
        if len(all_small_vols) > len(out_frac_small[ssp]):
            out_frac_small[ssp] = out_frac_small[ssp].reindex(
                all_small_vols.index)
            out_abs_small[ssp] = out_abs_small[ssp].reindex(
                all_small_vols.index)

        # use longest time index for regional volume dictionaries
        if len(all_large_vols) > len(out_frac_large[ssp]):
            out_frac_large[ssp] = out_frac_large[ssp].reindex(
                all_large_vols.index)
            out_abs_large[ssp] = out_abs_large[ssp].reindex(
                all_large_vols.index)

        # add to regional volume dictionaries
        out_frac_small[ssp][gcm] = all_small_vols / all_small_vols.loc[
            ref_year]
        out_frac_large[ssp][gcm] = all_large_vols / all_large_vols.loc[
            ref_year]
        out_abs_small[ssp][gcm] = all_small_vols
        out_abs_large[ssp][gcm] = all_large_vols

    # store as *.csv files
    for k, v in out_frac_small.items():
        v.to_csv(odir_frac + k + '_small.csv')
    for k, v in out_frac_large.items():
        v.to_csv(odir_frac + k + '_large.csv')
    for k, v in out_abs_small.items():
        v.to_csv(odir_abs + k + '_small.csv')
    for k, v in out_abs_large.items():
        v.to_csv(odir_abs + k + '_large.csv')


def merge_flowline():
    # specify region
    rgi_reg = int(os.environ['RGI_REG'])
    ref_year = 2020

    # define input directory
    idir = '/home/www/fmaussion/run_moritz_cmip/cmip6_output/'

    # get list of failed glaciers
    failed = pd.read_csv('/home/users/moberrauch/cmip6_run/'
                         'failed_glaciers_fl.csv',
                         index_col=0)
    failed = failed.values.reshape(len(failed))

    # specify path to output *.nc files for historical runs
    f_past_tpl = '/home/www/oggm/gdirs/oggm_v1.4/exps/thesis_vas/' \
                 'RGI62/b_080/L5/summary/historical_run_output_extended_{}.nc'

    # specify dummy path to glacier statistics
    stas_dirpath = '/home/www/oggm/gdirs/oggm_v1.4/exps/thesis_vas/' \
                   'RGI62/b_080/L5/summary/glacier_statistics_{}.csv'

    glacier_stats = pd.read_csv(stas_dirpath.format(f'{rgi_reg:02d}'),
                                index_col=0, low_memory=False)

    # some logging info
    print(f'RGI region: {rgi_reg}', flush=True)

    # relative volume of given region, in realtion to reference year
    out_frac_small = dict()
    out_frac_large = dict()
    # absolute volume of given region
    out_abs_small = dict()
    out_abs_large = dict()

    # get RGI IDs for large glaciers
    quantile = 0.95
    area_threshold = glacier_stats.rgi_area_km2.quantile(quantile)
    small_ids = glacier_stats[
        (glacier_stats.rgi_area_km2 <= area_threshold)].index.values
    large_ids = glacier_stats[
        (glacier_stats.rgi_area_km2 > area_threshold)].index.values

    # make sure the RGI ID is a string
    rgi_reg = '{:02d}'.format(rgi_reg)
    # get RGI IDs for failed glaciers (see section above)
    sel_ids = [r for r in failed if rgi_reg + '.' in r]

    sel_large_ids = np.unique(np.concatenate([sel_ids, small_ids]))
    sel_small_ids = np.unique(np.concatenate([sel_ids, large_ids]))

    # get path to historic run for given RGI region
    f_past = f_past_tpl.format(rgi_reg)

    # define paths and create output directories
    odir_frac = '/home/users/moberrauch/run_output/' \
                'merge_cru_cmip_fl/frac_output/' + rgi_reg + '/'
    odir_abs = '/home/users/moberrauch/run_output/' \
               'merge_cru_cmip_fl/abs_output/' + rgi_reg + '/'
    mkdir(odir_frac)
    mkdir(odir_abs)

    # get path to RGI region and the *.nc files
    # for all the GCMs and SSPs
    rdir = idir + 'RGI' + rgi_reg + '/'
    all_ncs = glob.glob(rdir + '*/*.nc')

    # select certain ssps
    ncs = list()
    for ssp_ in ['ssp126', 'ssp245', 'ssp370', 'ssp585']:
        ncs.extend([nc for nc in all_ncs if ssp_ in nc])

    # iterate over all *.nc files, i.e. all combinations of region/GCM/SSP
    for nc in progressbar.progressbar(sorted(ncs)):

        # get SSP and GCM description
        ssp = os.path.basename(nc).split('_')[1].replace('.nc', '')
        gcm = os.path.basename(nc).split('_')[0]
        # define key for dataset
        key_both = ssp + '_' + gcm
        # create new dataframe for give SSP, if it doesn't already exist
        # and add it to nested dictionary for absolute and relative volume
        if ssp not in out_frac_small:
            out_frac_small[ssp] = pd.DataFrame()
            out_frac_large[ssp] = pd.DataFrame()
            out_abs_small[ssp] = pd.DataFrame()
            out_abs_large[ssp] = pd.DataFrame()

        # merge historical run and forward run
        ds_small = extended_fl(f_past, nc, sel_ids=sel_small_ids)
        ds_large = extended_fl(f_past, nc, sel_ids=sel_large_ids)

        # sum volume over all glaciers and convert into pd.Series
        all_small_vols = ds_small.volume.sum(dim='rgi_id').to_series()
        all_large_vols = ds_large.volume.sum(dim='rgi_id').to_series()

        # use longest time index for regional volume dictionaries
        if len(all_small_vols) > len(out_frac_small[ssp]):
            out_frac_small[ssp] = out_frac_small[ssp].reindex(
                all_small_vols.index)
            out_abs_small[ssp] = out_abs_small[ssp].reindex(
                all_small_vols.index)

        # use longest time index for regional volume dictionaries
        if len(all_large_vols) > len(out_frac_large[ssp]):
            out_frac_large[ssp] = out_frac_large[ssp].reindex(
                all_large_vols.index)
            out_abs_large[ssp] = out_abs_large[ssp].reindex(
                all_large_vols.index)

        # add to regional volume dictionaries
        out_frac_small[ssp][gcm] = all_small_vols / all_small_vols.loc[
            ref_year]
        out_frac_large[ssp][gcm] = all_large_vols / all_large_vols.loc[
            ref_year]
        out_abs_small[ssp][gcm] = all_small_vols
        out_abs_large[ssp][gcm] = all_large_vols

    # store as *.csv files
    for k, v in out_frac_small.items():
        v.to_csv(odir_frac + k + '_small.csv')
    for k, v in out_frac_large.items():
        v.to_csv(odir_frac + k + '_large.csv')
    for k, v in out_abs_small.items():
        v.to_csv(odir_abs + k + '_small.csv')
    for k, v in out_abs_large.items():
        v.to_csv(odir_abs + k + '_large.csv')


if __name__ == '__main__':
    merge_vas()
    merge_flowline()

