import os
import logging
import sys

# Libs
import xarray as xr
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt

# Locals
import oggm.cfg as cfg
from oggm import utils, workflow, tasks
from oggm.core import gcm_climate

if __name__ == '__main__':

    # Initialize OGGM and set up the default run parameters
    cfg.initialize(logging_level='ERROR')
    rgi_version = '62'

    cfg.PARAMS['border'] = 80

    # Local working directory (where OGGM will write its output)
    WORKING_DIR = os.environ.get('OGGM_WORKDIR', '')
    if not WORKING_DIR:
        raise RuntimeError('Need a working dir')
    utils.mkdir(WORKING_DIR)
    cfg.PATHS['working_dir'] = WORKING_DIR

    OUTPUT_DIR = os.environ.get('OGGM_OUTDIR', '')
    if not OUTPUT_DIR:
        raise RuntimeError('Need an output dir')
    utils.mkdir(OUTPUT_DIR)

    cfg.PARAMS['prcp_scaling_factor'] = 1.8
    cfg.PARAMS['continue_on_error'] = True

    # Init
    workflow.init_mp_pool(True)

    rgi_reg = os.environ.get('OGGM_RGI_REG', '')
    if rgi_reg not in ['{:02d}'.format(r) for r in range(1, 20)]:
        raise RuntimeError('Need an RGI Region')

    # Module logger
    log = logging.getLogger(__name__)
    log.workflow('Starting run for RGI reg {}'.format(rgi_reg))

    # RGI glaciers
    rgi_ids = gpd.read_file(
        utils.get_rgi_region_file(rgi_reg, version=rgi_version))

    # For greenland we omit connectivity level 2
    if rgi_reg == '05':
        rgi_ids = rgi_ids.loc[rgi_ids['Connect'] != 2]

    # Go - get the pre-processed glacier directories
    base_url = 'https://cluster.klima.uni-bremen.de/~fmaussion/gdirs/final_prepro_cmip6/era5_eb'
    gdirs = workflow.init_glacier_directories(rgi_ids, from_prepro_level=5,
                                              prepro_base_url=base_url,
                                              prepro_rgi_version=rgi_version)

    gcms = pd.read_csv('/home/www/oggm/cmip6/all_gcm_list.csv', index_col=0)

    n_gcms = len(sys.argv) - 1

    for gcm in sys.argv[1:]:
        df1 = gcms.loc[gcms.gcm == gcm]
        for ssp in df1.ssp.unique():
            df2 = df1.loc[df1.ssp == ssp]
            assert len(df2) == 2
            ft = df2.loc[df2['var'] == 'tas'].iloc[0]
            fp = df2.loc[df2['var'] == 'pr'].iloc[0].path
            rid = ft.fname.replace('_r1i1p1f1_tas.nc', '')
            ft = ft.path

            log.workflow('Starting run for {}'.format(rid))

            workflow.execute_entity_task(gcm_climate.process_cmip5_data, gdirs,
                                         filesuffix='_' + rid,
                                         # recognize the climate file for later
                                         fpath_temp=ft,
                                         # temperature projections
                                         fpath_precip=fp,  # precip projections
                                         year_range=('1981', '2018'),
                                         );
            workflow.execute_entity_task(tasks.run_from_climate_data, gdirs,
                                         climate_filename='gcm_data',
                                         # use gcm_data, not climate_historical
                                         climate_input_filesuffix='_' + rid,
                                         # use a different scenario
                                         init_model_filesuffix='_historical',
                                         # this is important! Start from 2019 glacier
                                         output_filesuffix=rid,
                                         # recognize the run for later
                                         return_value=False,
                                         );
            gcm_dir = os.path.join(OUTPUT_DIR, 'RGI' + rgi_reg, gcm)
            utils.mkdir(gcm_dir)
            utils.compile_run_output(gdirs, input_filesuffix=rid,
                                     path=os.path.join(gcm_dir, rid + '.nc'))

    log.workflow('OGGM Done')

