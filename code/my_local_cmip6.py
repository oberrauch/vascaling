import os
import logging
import sys

# Libs
import xarray as xr
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt

# Locals
from oggm import cfg, utils, workflow, tasks
from oggm.core import gcm_climate

import oggm_vas as vascaling


def run_cmip():
    # Initialize OGGM and set up the default run parameters
    vascaling.initialize(logging_level='WORKFLOW')
    rgi_version = '62'
    cfg.PARAMS['border'] = 80

    # paths to working directory and output directory
    local_wdir = '/Users/oberrauch/work/master/working_directories/cmip6/'
    local_outdir = '/Users/oberrauch/work/master/data/cmip6_data/'
    wdir = os.environ.get('WORKDIR', local_wdir)
    cfg.PATHS['working_dir'] = wdir
    outdir = os.environ.get('OUTDIR', local_outdir)

    # define the baseline climate CRU or HISTALP
    cfg.PARAMS['baseline_climate'] = 'CRU'
    # set the mb hyper parameters accordingly
    cfg.PARAMS['prcp_scaling_factor'] = 3
    cfg.PARAMS['temp_melt'] = 0
    cfg.PARAMS['temp_all_solid'] = 4
    cfg.PARAMS['run_mb_calibration'] = False
    # set minimum ice thickness to include in glacier length computation
    # this reduces weird spikes in length records
    cfg.PARAMS['min_ice_thick_for_length'] = 0.1

    # the bias is defined to be zero during the calibration process,
    # which is why we don't use it here to reproduce the results
    cfg.PARAMS['use_bias_for_run'] = True

    # read RGI entry for the glaciers as DataFrame
    # containing the outline area as shapefile
    # RGI glaciers
    rgi_reg = '11'
    rgi_ids = gpd.read_file(
        utils.get_rgi_region_file(rgi_reg, version=rgi_version))

    rgi_ids = ['RGI60-11.00897']

    # get and set path to intersect shapefile
    intersects_db = utils.get_rgi_intersects_region_file(region=rgi_reg)
    cfg.set_intersects_db(intersects_db)

    # operational run, all glaciers should run
    cfg.PARAMS['continue_on_error'] = False


    # Module logger
    log = logging.getLogger(__name__)
    log.workflow('Starting run for RGI reg {}'.format(rgi_reg))

    # Go - get the pre-processed glacier directories
    base_url = 'https://cluster.klima.uni-bremen.de/~moberrauch/prepro/'
    gdirs = workflow.init_glacier_directories(rgi_ids, from_prepro_level=3,
                                              prepro_base_url=base_url,
                                              prepro_rgi_version=rgi_version)

    # run vascaling climate tasks
    # workflow.execute_entity_task(vascaling.local_t_star, gdirs)
    # adjust mass balance residual with geodetic observations
    # vascaling.match_regional_geodetic_mb(gdirs=gdirs, rgi_reg=rgi_reg)
    # prepare historic "spinup"
    # workflow.execute_entity_task(vascaling.run_from_climate_data, gdirs,
    #                              ys=2003, ye=2020,
    #                              output_filesuffix='_historical')
    # read gcm list
    gcms = pd.read_csv('https://cluster.klima.uni-bremen.de/'
                       + '~oggm/cmip6/all_gcm_list.csv', index_col=0)

    utils.gdir_to_tar

    # iterate over all specified gcms
    argv = ['CESM2']
    # iterate over all specified gcms
    for gcm in argv:
        # iterate over all SSPs (Shared Socioeconomic Pathways)
        df1 = gcms.loc[gcms.gcm == gcm]
        for ssp in df1.ssp.unique():
            df2 = df1.loc[df1.ssp == ssp]
            assert len(df2) == 2
            # get temperature projections
            ft = df2.loc[df2['var'] == 'tas'].iloc[0]
            # get precipitation projections
            fp = df2.loc[df2['var'] == 'pr'].iloc[0].path
            rid = ft.fname.replace('_r1i1p1f1_tas.nc', '')
            ft = ft.path

            # adjust paths for access from local
            ft = '/Users/oberrauch/Downloads/CESM2_ssp585_r1i1p1f1_tas.nc'
            fp = '/Users/oberrauch/Downloads/CESM2_ssp585_r1i1p1f1_pr.nc'

            log.workflow('Starting run for {}'.format(rid))

            workflow.execute_entity_task(gcm_climate.process_cmip_data, gdirs,
                                         filesuffix='_' + rid,
                                         # recognize the climate file for later
                                         fpath_temp=ft,
                                         # temperature projections
                                         fpath_precip=fp,  # precip projections
                                         year_range=('1981', '2020'))

            workflow.execute_entity_task(vascaling.run_from_climate_data,
                                         gdirs, climate_filename='gcm_data',
                                         climate_input_filesuffix='_' + rid,
                                         init_model_filesuffix='_historical',
                                         output_filesuffix=rid,
                                         return_value=False)
            gcm_dir = os.path.join(outdir, 'RGI' + rgi_reg, gcm)
            utils.mkdir(gcm_dir)
            utils.compile_run_output(gdirs, input_filesuffix=rid,
                                     path=os.path.join(gcm_dir, rid + '.nc'))


    log.workflow('OGGM Done')


if __name__ == '__main__':
    run_cmip()
