""" The `equilibrium_run.py` script combines the run output of the different
models into one :py:class:`xarray.Dataset`. This happens at the end of the
script, after all runs are completed. Hence, if one run fails nothing happens.
Therefore, I decided to skip this part on the cluster run and store the files
separately. Now I have to combine them...
"""

import os
import shutil
import numpy as np
import pandas as pd
import xarray as xr

import logging
logging.basicConfig(level=logging.INFO)

log = logging.getLogger('combine output')

# list all files in directory
in_dir = '/home/users/moberrauch/run_output/'
files = os.listdir(in_dir)
log.info('Files in the input folder: {}'.format(', '.join(files)))

# output directory
reset = False
out_dir = '/home/users/moberrauch/run_output/'
if not os.path.exists(out_dir) or reset:
    log.info('Creating output directory')
    try:
        shutil.rmtree(out_dir)
    except:
        raise
    os.mkdir(out_dir)

# define models and mb model suffixes
models = ['fl', 'vas']
mb_models = ['constant', 'random']

# combine the mass balance datasets
# ---------------------------------

# open datasets
ds = list()
for model in models:
    log.info('Reading {} mass balance files'.format(model))
    ds_ = xr.open_dataset(os.path.join(in_dir,
                                       'mb_output_{:s}.nc').format(model))
    ds.append(ds_)

# concat datasets by evolution model
log.info('Combining mass balance datasets')
ds = xr.concat(ds, pd.Index(models, name='model'))
# store to file
out_path = os.path.join(out_dir, 'mb_output.nc')
log.info('Store combined mass balance dataset to {}'.format(out_path))
ds.to_netcdf(out_path)

# close datasets
ds_.close()
ds.close()

# combine the model output datasets
# ---------------------------------

# open datasets
ds = list()
for mb_model in mb_models:
    ds_ = list()
    for model in models:
        log.info('Reading {} mb, {} run output file'.format(mb_model, model))
        f_name = 'run_output_{:s}_{:s}.nc'.format(mb_model, model)
        ds_.append(xr.open_dataset(os.path.join(in_dir, f_name)))

    # concat datasets by model
    log.info('Combining {} mb model by evolution model'. format(mb_model))
    ds.append(xr.concat(ds_, 'model'))

# concat datasets by mass balance model
log.info('Combining by mass balance')
ds = xr.concat(ds, 'mb_model')

# store to file
out_path = os.path.join(out_dir, 'eq_runs.nc')
log.info('Store combined run output dataset to {}'.format(out_path))
ds.to_netcdf(out_path)

