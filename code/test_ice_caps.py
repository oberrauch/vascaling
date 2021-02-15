# import internal and externals libraries
import os
import shutil
import numpy as np
import pandas as pd
import xarray as xr

import logging

log = logging.getLogger('test-ice-caps')

# import the needed OGGM modules
from oggm import cfg, utils, workflow
from oggm.core import gis, climate, flowline
import oggm_vas as vascaling

log.info('Starting run')

# specify glaciers by RGI IDs (INPUT)
rgi_ids = ['RGI60-15.00105', 'RGI60-15.00107', 'RGI60-15.00194',
           'RGI60-15.01653', 'RGI60-15.01654', 'RGI60-15.01658',
           'RGI60-15.02072', 'RGI60-15.02380', 'RGI60-15.02710',
           'RGI60-15.02817', 'RGI60-15.02841', 'RGI60-15.02842',
           'RGI60-15.02864', 'RGI60-15.02909', 'RGI60-15.02956',
           # 'RGI60-15.03030', 'RGI60-15.03053', 'RGI60-15.03056',
           # 'RGI60-15.03109', 'RGI60-15.03211', 'RGI60-15.03273',
           # 'RGI60-15.03276', 'RGI60-15.03337', 'RGI60-15.03404',
           # 'RGI60-15.03586', 'RGI60-15.03591', 'RGI60-15.03592',
           # 'RGI60-15.03613', 'RGI60-15.03722', 'RGI60-15.03729',
           # 'RGI60-15.03967', 'RGI60-15.04154', 'RGI60-15.04163',
           # 'RGI60-15.04278', 'RGI60-15.04338', 'RGI60-15.04584',
           # 'RGI60-15.04590', 'RGI60-15.04629', 'RGI60-15.04630',
           # 'RGI60-15.04659', 'RGI60-15.04760', 'RGI60-15.04855',
           # 'RGI60-15.04863', 'RGI60-15.04866', 'RGI60-15.04996',
           # 'RGI60-15.05121', 'RGI60-15.05515', 'RGI60-15.05538',
           # 'RGI60-15.07062', 'RGI60-15.07265', 'RGI60-15.07601',
           # 'RGI60-15.07736', 'RGI60-15.07738', 'RGI60-15.07739',
           # 'RGI60-15.07830', 'RGI60-15.07836', 'RGI60-15.07868',
           # 'RGI60-15.07869', 'RGI60-15.07872', 'RGI60-15.07877',
           # 'RGI60-15.07880', 'RGI60-15.08168', 'RGI60-15.08195',
           # 'RGI60-15.08200', 'RGI60-15.08282', 'RGI60-15.09121',
           # 'RGI60-15.09181', 'RGI60-15.09227', 'RGI60-15.09260',
           # 'RGI60-15.09269', 'RGI60-15.09271', 'RGI60-15.09272',
           # 'RGI60-15.09283', 'RGI60-15.09289', 'RGI60-15.09331',
           # 'RGI60-15.09339', 'RGI60-15.09341', 'RGI60-15.09354',
           # 'RGI60-15.09361', 'RGI60-15.09377', 'RGI60-15.09380',
           # 'RGI60-15.09383', 'RGI60-15.09438', 'RGI60-15.09444',
           # 'RGI60-15.09453', 'RGI60-15.09456', 'RGI60-15.09487',
           # 'RGI60-15.09758', 'RGI60-15.09795', 'RGI60-15.09799',
           # 'RGI60-15.09802', 'RGI60-15.09808', 'RGI60-15.09835',
           # 'RGI60-15.09847', 'RGI60-15.09950', 'RGI60-15.10060',
           # 'RGI60-15.10070', 'RGI60-15.10076', 'RGI60-15.10115',
           # 'RGI60-15.10190', 'RGI60-15.10192', 'RGI60-15.10222',
           # 'RGI60-15.10251', 'RGI60-15.10275', 'RGI60-15.10286',
           # 'RGI60-15.10288', 'RGI60-15.10303', 'RGI60-15.10313',
           # 'RGI60-15.10314', 'RGI60-15.10399', 'RGI60-15.10428',
           # 'RGI60-15.10444', 'RGI60-15.10564', 'RGI60-15.10567',
           # 'RGI60-15.10569', 'RGI60-15.10577', 'RGI60-15.10727',
           # 'RGI60-15.10737', 'RGI60-15.10741', 'RGI60-15.10764',
           # 'RGI60-15.10765', 'RGI60-15.10821', 'RGI60-15.10826',
           # 'RGI60-15.10847', 'RGI60-15.10879', 'RGI60-15.10883',
           # 'RGI60-15.10885', 'RGI60-15.10886', 'RGI60-15.10887',
           # 'RGI60-15.10890', 'RGI60-15.10895', 'RGI60-15.10923',
           # 'RGI60-15.10928', 'RGI60-15.10929', 'RGI60-15.10932',
           # 'RGI60-15.10935', 'RGI60-15.10989', 'RGI60-15.11015',
           # 'RGI60-15.11028', 'RGI60-15.11145', 'RGI60-15.11154',
           # 'RGI60-15.11164', 'RGI60-15.11560', 'RGI60-15.11562',
           # 'RGI60-15.11563', 'RGI60-15.11604', 'RGI60-15.11619',
           # 'RGI60-15.11620', 'RGI60-15.11622', 'RGI60-15.11696',
           # 'RGI60-15.11701', 'RGI60-15.11709', 'RGI60-15.11758',
           # 'RGI60-15.11759', 'RGI60-15.11763', 'RGI60-15.11770',
           # 'RGI60-15.11771', 'RGI60-15.11782', 'RGI60-15.11785',
           # 'RGI60-15.11787', 'RGI60-15.11788', 'RGI60-15.11789',
           # 'RGI60-15.11790', 'RGI60-15.11795', 'RGI60-15.11798',
           # 'RGI60-15.11808', 'RGI60-15.11812', 'RGI60-15.11855',
           # 'RGI60-15.11867', 'RGI60-15.11879', 'RGI60-15.11888',
           # 'RGI60-15.11909', 'RGI60-15.11926', 'RGI60-15.11930',
           # 'RGI60-15.11957', 'RGI60-15.11975', 'RGI60-15.11988',
           # 'RGI60-15.11989', 'RGI60-15.12033', 'RGI60-15.12040',
           # 'RGI60-15.12041', 'RGI60-15.12067', 'RGI60-15.12084',
           # 'RGI60-15.12090', 'RGI60-15.12091', 'RGI60-15.12092',
           # 'RGI60-15.12099', 'RGI60-15.12100', 'RGI60-15.12104',
           # 'RGI60-15.12109', 'RGI60-15.12111', 'RGI60-15.12112',
           # 'RGI60-15.12114', 'RGI60-15.12116', 'RGI60-15.12121',
           # 'RGI60-15.12122', 'RGI60-15.12128', 'RGI60-15.12151',
           # 'RGI60-15.12152', 'RGI60-15.12153', 'RGI60-15.12166',
           # 'RGI60-15.12183', 'RGI60-15.12185', 'RGI60-15.12210',
           # 'RGI60-15.12212', 'RGI60-15.12215', 'RGI60-15.12217',
           # 'RGI60-15.12242', 'RGI60-15.12243', 'RGI60-15.12246',
           # 'RGI60-15.12247', 'RGI60-15.12248', 'RGI60-15.12249',
           # 'RGI60-15.12258', 'RGI60-15.12263', 'RGI60-15.12266',
           # 'RGI60-15.12284', 'RGI60-15.12289', 'RGI60-15.12291',
           # 'RGI60-15.12298', 'RGI60-15.12304', 'RGI60-15.12314',
           # 'RGI60-15.12339', 'RGI60-15.12365', 'RGI60-15.12411',
           # 'RGI60-15.12412', 'RGI60-15.12415', 'RGI60-15.12434',
           # 'RGI60-15.12435', 'RGI60-15.12438', 'RGI60-15.12500',
           # 'RGI60-15.12502', 'RGI60-15.12505', 'RGI60-15.12506',
           # 'RGI60-15.12520', 'RGI60-15.12544', 'RGI60-15.12547',
           # 'RGI60-15.12549', 'RGI60-15.12551', 'RGI60-15.12552',
           # 'RGI60-15.12553', 'RGI60-15.12566', 'RGI60-15.12586',
           # 'RGI60-15.12594', 'RGI60-15.12627', 'RGI60-15.12638',
           # 'RGI60-15.12644', 'RGI60-15.12654', 'RGI60-15.12667',
           # 'RGI60-15.12686', 'RGI60-15.12693', 'RGI60-15.12723',
           # 'RGI60-15.12726', 'RGI60-15.12734', 'RGI60-15.12737',
           # 'RGI60-15.12740', 'RGI60-15.12754', 'RGI60-15.12775',
           # 'RGI60-15.12792', 'RGI60-15.12793', 'RGI60-15.12809',
           # 'RGI60-15.12829', 'RGI60-15.12833', 'RGI60-15.12880',
           # 'RGI60-15.12881', 'RGI60-15.12889', 'RGI60-15.12909',
           # 'RGI60-15.12921', 'RGI60-15.12922', 'RGI60-15.12945',
           # 'RGI60-15.12958', 'RGI60-15.12959', 'RGI60-15.12960',
           # 'RGI60-15.12975', 'RGI60-15.12979', 'RGI60-15.13009',
           # 'RGI60-15.13013', 'RGI60-15.13033', 'RGI60-15.13035',
           # 'RGI60-15.13063', 'RGI60-15.13066', 'RGI60-15.13067',
           # 'RGI60-15.13070', 'RGI60-15.13071', 'RGI60-15.13082',
           # 'RGI60-15.13083', 'RGI60-15.13088', 'RGI60-15.13091',
           # 'RGI60-15.13092'
           ]

# compute RGI region and version from RGI IDs
# assuming all they are all the same
rgi_region = (rgi_ids[0].split('-')[-1]).split('.')[0]
rgi_version = '62'

# load default parameter file
vascaling.initialize()

# # get path to directories on the CLUSTER - comment/uncomment as necessary
# OUTPUT_DIR = os.environ['OUTDIR']
# WORKING_DIR = os.environ['WORKDIR']

# get LOCAL environmental variables for working and output directories
WORKING_DIR = '/Users/oberrauch/work/master/working_directories/ice-caps-new/'
# create working directory
utils.mkdir(WORKING_DIR)

# specify path to working directory
cfg.PATHS['working_dir'] = WORKING_DIR
# define the baseline climate CRU or HISTALP
cfg.PARAMS['baseline_climate'] = 'CRU'
# set the mb hyper parameters accordingly
cfg.PARAMS['prcp_scaling_factor'] = 3
cfg.PARAMS['temp_melt'] = 0
cfg.PARAMS['temp_all_solid'] = 4
cfg.PARAMS['run_mb_calibration'] = False
# set map size
cfg.PARAMS['border'] = 80
# set minimum ice thickness to include in glacier length computation
# this reduces weird spikes in length records
cfg.PARAMS['min_ice_thick_for_length'] = 0.1

# the bias is defined to be zero during the calibration process,
# which is why we don't use it here to reproduce the results
cfg.PARAMS['use_bias_for_run'] = True

# get and set path to intersect shapefile
intersects_db = utils.get_rgi_intersects_region_file(region=rgi_region)
cfg.set_intersects_db(intersects_db)

# operational run, all glaciers should run
cfg.PARAMS['continue_on_error'] = True

# Go - get the pre-processed glacier directories
base_url = 'https://cluster.klima.uni-bremen.de/~oggm/gdirs/oggm_v1.4/' \
           'L3-L5_files/RGIV62_fleb_qc3_CRU_pcp2.5'
gdirs = workflow.init_glacier_directories(rgi_ids, from_prepro_level=3,
                                          prepro_base_url=base_url,
                                          prepro_rgi_version=rgi_version)

# compute local t* and the corresponding mu*
workflow.execute_entity_task(vascaling.local_t_star, gdirs)

# run and compile output
workflow.execute_entity_task(vascaling.run_from_climate_data, gdirs)
utils.compile_run_output(gdirs)
