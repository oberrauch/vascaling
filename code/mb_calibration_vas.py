""" Run the mass balance calibration for the VAS model. """

# Python imports
import json
import os

# Locals
from oggm import cfg, utils, tasks, workflow
from oggm.core import gis
from oggm.workflow import execute_entity_task
import oggm_vas as vascaling


def mb_calibration(rgi_version, baseline):
    """ Run the mass balance calibration for the VAS model. RGI version and
    baseline cliamte must be given.

    :param rgi_version: int, RGI version
    :param baseline: str, baseline climate 'HISTALP' or 'CRU'
    """

    # initialize OGGM and set up the run parameters
    vascaling.initialize(logging_level='WORKFLOW')

    # LOCAL paths (where to write the OGGM run output)
    # dirname = 'VAS_ref_mb_{}_RGIV{}'.format(baseline, rgi_version)
    # wdir = utils.gettempdir(dirname, home=True, reset=True)
    # utils.mkdir(wdir, reset=True)
    # cfg.PATHS['working_dir'] = wdir

    # CLUSTER paths
    wdir = os.environ.get('WORKDIR', '')
    cfg.PATHS['working_dir'] = wdir

    # we are running the calibration ourselves
    cfg.PARAMS['run_mb_calibration'] = True
    # we are using which baseline data?
    cfg.PARAMS['baseline_climate'] = baseline
    # no need for intersects since this has an effect on the inversion only
    cfg.PARAMS['use_intersects'] = False
    # use multiprocessing?
    cfg.PARAMS['use_multiprocessing'] = True
    # set to True for operational runs
    cfg.PARAMS['continue_on_error'] = True
    # 10 is only for OGGM-VAS, OGGM needs 80 to run
    cfg.PARAMS['border'] = 80

    if baseline == 'HISTALP':
        # other params: see https://oggm.org/2018/08/10/histalp-parameters/
        # cfg.PARAMS['prcp_scaling_factor'] = 1.75
        # cfg.PARAMS['temp_melt'] = -1.75
        cfg.PARAMS['prcp_scaling_factor'] = 2.5
        cfg.PARAMS['temp_melt'] = -0.5
    elif baseline == 'CRU':
        # using the parameters from Marzeion et al. (2012)
        # cfg.PARAMS['prcp_scaling_factor'] = 2.5
        # cfg.PARAMS['temp_melt'] = 1
        # cfg.PARAMS['temp_all_solid'] = 3
        # using the parameters from Malles and Marzeion 2020
        cfg.PARAMS['prcp_scaling_factor'] = 3
        cfg.PARAMS['temp_melt'] = 0
        cfg.PARAMS['temp_all_solid'] = 4
        # cfg.PARAMS['prcp_gradient'] = 4

    # the next step is to get all the reference glaciers,
    # i.e. glaciers with mass balance measurements.

    # get the reference glacier ids (they are different for each RGI version)
    df, _ = utils.get_wgms_files()
    rids = df['RGI{}0_ID'.format(rgi_version[0])]

    # we can't do Antarctica
    rids = [rid for rid in rids if not ('-19.' in rid)]

    # For HISTALP only RGI reg 11.01 (ALPS)
    if baseline == 'HISTALP':
        rids = [rid for rid in rids if '-11' in rid]

    # make a new dataframe with those (this takes a while)
    print('Reading the RGI shapefiles...')
    rgidf = utils.get_rgi_glacier_entities(rids, version=rgi_version)
    print('For RGIV{} we have {} candidate reference '
          'glaciers.'.format(rgi_version, len(rgidf)))

    # initialize the glacier regions
    base_url = 'https://cluster.klima.uni-bremen.de/~oggm/gdirs/oggm_v1.4/' \
               'L3-L5_files/RGIV62_fleb_qc3_CRU_pcp2.5'
    # Go - get the pre-processed glacier directories
    gdirs = workflow.init_glacier_directories(rids, from_prepro_level=3,
                                              prepro_base_url=base_url,
                                              prepro_rgi_version=rgi_version)

    # Some glaciers in RGI Region 11 are not inside the HISTALP domain
    if baseline == 'HISTALP':
        gdirs = [gdir for gdir in gdirs if gdir.rgi_subregion == '11-01']

    # get reference glaciers with mass balance measurements
    gdirs = utils.get_ref_mb_glaciers(gdirs)

    # keep only these glaciers
    rgidf = rgidf.loc[rgidf.RGIId.isin([g.rgi_id for g in gdirs])]

    # save to file
    rgidf.to_file(os.path.join(wdir, 'mb_ref_glaciers.shp'))
    print('For RGIV{} and {} we have {} reference glaciers'.format(rgi_version,
                                                                   baseline,
                                                                   len(rgidf)))

    # sort for more efficient parallel computing
    rgidf = rgidf.sort_values('Area', ascending=False)

    # newly initialize glacier directories
    gdirs = workflow.init_glacier_directories(rgidf, reset=False, force=True)
    workflow.execute_entity_task(gis.define_glacier_region, gdirs)
    workflow.execute_entity_task(gis.glacier_masks, gdirs)

    # run climate tasks
    vascaling.compute_ref_t_stars(gdirs)
    execute_entity_task(vascaling.local_t_star, gdirs)

    # we store the associated params
    mb_calib = gdirs[0].read_pickle('climate_info')['mb_calib_params']
    with open(os.path.join(wdir, 'mb_calib_params.json'), 'w') as fp:
        json.dump(mb_calib, fp)


if __name__ == '__main__':
    """ Run the VAS model mass balance calibration for the given RGI version
    and baseline climate file. `ref_tstars.csv` is store in
    ~/tmp/OGGM/VAS_ref_mb_xx_RGIVxx/ (locally) directory or
    ~/run_output/mb_calib_cru/ directory (on the cluster).
    """
    # specify RGI Version and baseline climate
    rgi_version = '62'
    baseline = 'CRU'

    # run mass balance calibration
    mb_calibration(rgi_version, baseline)
