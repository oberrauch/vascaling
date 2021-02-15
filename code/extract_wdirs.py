#!/usr/bin/env python
# coding: utf-8

# Extract glacier directories
import os
import geopandas as gpd
from oggm import cfg, utils, workflow

# define paths
prepro_dir = '/home/users/moberrauch/run_output/vas_prepro'
working_dir = '/home/users/moberrauch/wdirs/historical'

# OGGM initialization
cfg.initialize()
cfg.PATHS['working_dir'] = working_dir


# specify RGI version and regions
rgi_version = '62'
rgi_regions = [11, 13, 14, 15]

# get RGI IDs
for rgi_region in rgi_regions:
    fpath = utils.get_rgi_region_file(rgi_region, rgi_version)
    if rgi_region == rgi_regions[0]:
        rgi_ids = gpd.read_file(fpath)
    else:
        rgi_ids = rgi_ids.append(gpd.read_file(fpath))

# extract working directories from *.tar files
workflow.init_glacier_directories(rgi_ids, from_tar=prepro_dir)

