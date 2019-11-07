import os
import xarray as xr

def combine_data_sets(wdir_template):
    """ Combines run output datasets for different (vas vs. flowline) models.

    Parameters
    ----------
    wdir_template: string
        Path to the working directory, where the output files are stored.
        Should contain one format placeholder {}

    Returns
    -------

    """
    # use the vas model working directory first
    wdir = wdir_template.format('vas')
    
    # vas, random: specify file path, and load dataset
    fpath = os.path.join(wdir, 'run_output_random_vas.nc')
    ds_vas_rand = xr.open_dataset(fpath)
    # vas, random, norm: specify file path, and load dataset
    fpath = os.path.join(wdir, 'normalized_output_random_vas.nc')
    ds_vas_norm_rand = xr.open_dataset(fpath)
    # vas, constans: specify file path, and load dataset
    fpath = os.path.join(wdir, 'run_output_constant_vas.nc')
    ds_vas_const = xr.open_dataset(fpath)
    # vas, constant, norm: specify file path, and load dataset
    fpath = os.path.join(wdir, 'normalized_output_constant_vas.nc')
    ds_vas_norm_const = xr.open_dataset(fpath)

    # change to the flowline model working directory
    wdir = wdir_template.format('fl')

    # fl, random: specify file path, and load dataset
    fpath = os.path.join(wdir, 'run_output_random_fl.nc')
    ds_fl_rand = xr.open_dataset(fpath)
    # fl, random, norm: specify file path, and load dataset
    fpath = os.path.join(wdir, 'normalized_output_random_fl.nc')
    ds_fl_norm_rand = xr.open_dataset(fpath)
    # fl, constans: specify file path, and load dataset
    fpath = os.path.join(wdir, 'run_output_constant_fl.nc')
    ds_fl_const = xr.open_dataset(fpath)
    # fl, constant, norm: specify file path, and load dataset
    fpath = os.path.join(wdir, 'normalized_output_constant_fl.nc')
    ds_fl_norm_const = xr.open_dataset(fpath)

    # concat output datasets by model types
    ds_rand = xr.concat([ds_vas_rand, ds_fl_rand], 'model')
    ds_norm_rand = xr.concat([ds_vas_norm_rand, ds_fl_norm_rand], 'model')
    ds_const = xr.concat([ds_vas_const, ds_fl_const], 'model')
    ds_norm_const = xr.concat([ds_vas_norm_const, ds_fl_norm_const], 'model')


if __name__ == '__main__':
    # specify path
    wdir = '/Users/oberrauch/work/master/working_directories/' + \
           'equilibrium_{}_wdir/'
    # run
    combine_data_sets(wdir)
