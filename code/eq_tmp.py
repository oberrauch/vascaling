from oggm import cfg, utils, workflow
import xarray as xr

# specify suffixes
suffixes = ['_normal', '_bias_p', '_bias_n']

rgi_ids = ['RGI60-11.00897']

# compute RGI region and version from RGI IDs
# assuming all they are all the same
rgi_region = (rgi_ids[0].split('-')[-1]).split('.')[0]
rgi_version = (rgi_ids[0].split('-')[0])[-2:]

# load default parameter file
cfg.initialize()

# create working directory
wdir = '/Users/oberrauch/work/master/working_directories/'
wdir += 'equilibrium_fl_wdir'

# set path to working directory
cfg.PATHS['working_dir'] = wdir

rgidf = utils.get_rgi_glacier_entities(rgi_ids)
gdirs = workflow.init_glacier_regions(rgidf)

utils.compile_run_output(gdirs, filesuffix='_normal')