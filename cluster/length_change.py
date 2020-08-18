import numpy as np
import pandas as pd
import xarray as xr

# open dataset
in_path = '/home/users/moberrauch/run_output/eq_runs.nc'
ds = xr.open_dataset(in_path)

# select one run from dataset
ds = ds.sel(model='vas',
            mb_model='random',
            normalized=int(False),
            temp_bias=0.5)
# drop mean and sum
ds = ds.drop_sel(rgi_id=['mean', 'sum'])
# select variable
ds = ds.length

# iterate over all glaciers by rgi_id
length_change = list()
for rgi_id in ds.rgi_id:
    length_0 = ds.sel(rgi_id=rgi_id).isel(time=0)
    length_end = np.mean(ds.sel(rgi_id=rgi_id).isel(time=slice(-1000, None)))
    length_change.append((length_end-length_0).values)

# create pandas Series and store to file
length_change = pd.Series(length_change, ds.rgi_id)
f_path = '/home/users/moberrauch/length_change.csv'
length_change.to_csv(f_path)