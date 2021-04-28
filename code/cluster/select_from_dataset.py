import numpy as np
import pandas as pd
import xarray as xr

# specify path and open dataset
f_path = '/home/users/moberrauch/run_output/mb_output.nc'
ds = xr.open_dataset(f_path)

# select by RGI ID
#rgi_sel = ['RGI60-11.00002', 'RGI60-11.00012', 'RGI60-11.00073',
#           'RGI60-11.00080', 'RGI60-11.00106', 'RGI60-11.00190',
#           'RGI60-11.00251', 'RGI60-11.00289', 'RGI60-11.00300',
#           'RGI60-11.00603', 'RGI60-11.00619', 'RGI60-11.00638',
#           'RGI60-11.00647', 'RGI60-11.00719', 'RGI60-11.00781',
#           'RGI60-11.00787', 'RGI60-11.00797', 'RGI60-11.00804',
#           'RGI60-11.00807', 'RGI60-11.00843', 'RGI60-11.00892',
#           'RGI60-11.00897', 'RGI60-11.00918', 'RGI60-11.00929',
#           'RGI60-11.01238', 'RGI60-11.01450', 'RGI60-11.01662',
#           'RGI60-11.01704', 'RGI60-11.01776', 'RGI60-11.01834',
#           'RGI60-11.01876', 'RGI60-11.01930', 'RGI60-11.01987',
#           'RGI60-11.02072', 'RGI60-11.02214', 'RGI60-11.02245',
#           'RGI60-11.02249', 'RGI60-11.02285', 'RGI60-11.02648',
#           'RGI60-11.02671', 'RGI60-11.02679', 'RGI60-11.02704',
#           'RGI60-11.02746', 'RGI60-11.02764', 'RGI60-11.02766',
#           'RGI60-11.02773', 'RGI60-11.02774', 'RGI60-11.03084',
#           'RGI60-11.03135', 'RGI60-11.03166', 'RGI60-11.03246',
#           'RGI60-11.03638', 'RGI60-11.03643', 'RGI60-11.03671',
#           'RGI60-11.03674', 'RGI60-11.03756', 'mean', 'sum']
rgi_sel = ['RGI60-11.01450', 'RGI60-11.03643', 'RGI60-11.00106', 'RGI60-11.01238',
           'RGI60-11.02766', 'RGI60-11.03638', 'RGI60-11.02773', 'RGI60-11.02704',
           'RGI60-11.00719', 'RGI60-11.02072', 'RGI60-11.00897']
ds_sel = ds.sel(rgi_id=rgi_sel)

# store to file
out_path = '/home/users/moberrauch/run_output/mb_output_sel.nc'
ds_sel.to_netcdf(out_path)
