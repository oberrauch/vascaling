import numpy as np
import pandas as pd
import xarray as xr

# specify path and open dataset
# f_path = '/home/users/moberrauch/run_output/'
f_path = '/Users/oberrauch/work/master/cluster/run_output/eq_runs.nc'
ds = xr.open_dataset(f_path)

# select by RGI ID
rgi_sel = ['mean', 'sum']
ds_sel = ds.sel(rgi_id=rgi_sel)

# store to file
# out_path = '/home/users/moberrauch/run_output/eq_run_sel.nc'
out_path = '/Users/oberrauch/work/master/cluster/run_output/eq_run_sel.nc'
ds_sel.to_netcdf(out_path)
