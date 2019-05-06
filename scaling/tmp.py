# Run and find start ares
# Set up an OGGM/VAS run from scratch and test the start area seeking tasks.
# The only thing to specify ist the RGI ID and the glacier's name, the rest
# should run without any adjustments...

## Import section
# import externals libs
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# import the needed OGGM modules
import oggm
from oggm import cfg, utils
from oggm.utils import get_rgi_glacier_entities
from oggm.core import gis, climate, centerlines
from oggm.core import vascaling

## decide wheter to show plots or not
show = False

## Initialization
# load parameter file
cfg.initialize()

## RGI entity
# choose glacier
name = 'Ob. Grindelwaldgletscher'
rgi_id = 'RGI60-11.01270'

# get/downlaod the rgi entity including the outline shapefile
rgi_df = get_rgi_glacier_entities([rgi_id])
# set name, since not delivered with RGI
if rgi_df.loc[int(rgi_id[-5:])-1, 'Name'] is None:
    rgi_df.loc[int(rgi_id[-5:])-1, 'Name'] = name

# select single entry
rgi_entity = rgi_df.iloc[0]
# visualize
rgi_df.plot()
plt.title(rgi_entity.Name)
if show:
    plt.show()


## Glacier Directory
# specify the working directory and define the glacier directory
cfg.PATHS['working_dir'] = './'
gdir = oggm.GlacierDirectory(rgi_entity, reset=True)

## DEM and GIS tasks
# get the path to the DEM file (will download if necessary)
dem = utils.get_topo_file(gdir.cenlon, gdir.cenlat)
print('DEM source: {}, path to DEM file: {}'.format(dem[1], dem[0][0]))

# set path in config file
cfg.PATHS['dem_file'] = dem[0][0]
cfg.PARAMS['border'] = 10
cfg.PARAMS['use_intersects'] = False

# run GIS tasks
gis.define_glacier_region(gdir, entity=rgi_entity)
gis.glacier_masks(gdir)

## Climate data
# set mb calibration parameters
cfg.PARAMS['baseline_climate'] = 'HISTALP'
cfg.PARAMS['prcp_scaling_factor'] = 1.75
cfg.PARAMS['temp_melt'] = -1.75
# process HistAlp climate data
climate.process_histalp_data(gdir)

## Centerlines
# run center line preprocessing tasks
centerlines.compute_centerlines(gdir)
centerlines.initialize_flowlines(gdir)
centerlines.compute_downstream_line(gdir)
centerlines.compute_downstream_bedshape(gdir)
centerlines.catchment_area(gdir)
centerlines.catchment_intersections(gdir)
centerlines.catchment_width_geom(gdir)
centerlines.catchment_width_correction(gdir)

## Mass balance model
#  compute local t* and the corresponding mu*
vascaling.local_t_star(gdir)
# see calibration results
print(gdir.read_json('vascaling_mustar'))

# create mass balance model
mb_mod = vascaling.VAScalingMassBalance(gdir)

# look at specific mass balance over climate data period
min_hgt, max_hgt = vascaling.get_min_max_elevation(gdir)
y0 = 1802
y1 = 2014
years = np.arange(y0, y1)
mb = list()
for y in years:
    mb.append(mb_mod.get_specific_mb(min_hgt, max_hgt, y))

# visualize
plt.plot(years, mb)
plt.axhline(0, c='k', ls=':', lw=0.8)
plt.title('Modeled mass balance - {}'.format('Ob. Grindelwald Gletscher'))
plt.ylabel('Mass balance [mm w.e. yr$^{-1}$]')
if show:
    plt.show()


## Find start area
# run scalar minimization
minimize_res = vascaling.find_start_area(gdir)
print(minimize_res)

# stop script if minimization was not successful
if minimize_res.status:
    sys.exit(minimize_res.status)

# instance glacier with today's values
model_ref = vascaling.VAScalingModel(year_0=gdir.rgi_date,
                                     area_m2_0=gdir.rgi_area_m2,
                                     min_hgt=min_hgt, max_hgt=max_hgt,
                                     mb_model=mb_mod)

# instance guessed starting areas
num = 15
area_guess = np.linspace(100, gdir.rgi_area_m2*2,  num, endpoint=True)
# create empty containers
iteration_list = list()
spec_mb_list = list()

# iterate over all starting areas
for area_ in area_guess:
    # instance iteration model
    model_guess = vascaling.VAScalingModel(year_0=gdir.rgi_date,
                                           area_m2_0=gdir.rgi_area_m2,
                                           min_hgt=min_hgt, max_hgt=max_hgt,
                                           mb_model=mb_mod)
    # set new starting values
    model_guess.create_start_glacier(area_)
    # run model and save years and area
    diag_ds = model_guess.run_until_and_store(year_end=model_ref.year)
    # create series and store in container
    iteration_list.append(diag_ds.area_m2.to_dataframe()['area_m2'])
    spec_mb_list.append(diag_ds.spec_mb.to_dataframe()['spec_mb'])
    
# create DataFrame
iteration_df = pd.DataFrame(iteration_list, index=['{:.2f}'.format(a/1e6) for a in area_guess])
iteration_df.index.name = 'Start Area [km$^2$]'

# set up model with resulted starting area
model = vascaling.VAScalingModel(year_0=model_ref.year_0,
                                 area_m2_0=model_ref.area_m2_0,
                                 min_hgt=model_ref.min_hgt_0,
                                 max_hgt=model_ref.max_hgt,
                                 mb_model=model_ref.mb_model)
model.create_start_glacier(minimize_res.x)

# run
diag_ds = model.run_until_and_store(year_end=model_ref.year)

plt.plot([0,1], [1,3])
plt.show()

# create figure and add axes
fig = plt.figure(figsize=[8, 6])
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
plt.plot([0,1], [1,1])
# plot
# ax.plot(diag_ds.time, diag_ds.area_m2, color='k', ls='--', lw=1.2, label='Best Guess')
# ax.axhline(model_ref.area_m2_0, c='k', ls=':', label='measured area in {}'.format(gdir.rgi_date))
# ax = iteration_df.T.plot(legend=False, figsize=[8,6], colormap='Spectral', ax=ax)
# # add legend
# handels, labels = ax.get_legend_handles_labels()
# labels[2:] = [r'{} km$^2$'.format(l) for l in labels[2:]]
# ax.legend(handels, labels, bbox_to_anchor=(1.05, 0.5), loc=6)
#
# # replot best guess estimate (in case it lies below another guess)
# ax.plot(diag_ds.time, diag_ds.area_m2, color='k', ls='--', lw=1.2, label='Best Guess')
#
# # labels, title
# ax.set_ylabel('Glacier area [m$^2$]')
# ax.set_title('Modelled glacier area - {}'.format(rgi_entity.Name))

#plt.gcf().savefig('Users/oberrauch/Desktop/start_area.png', bbox_inches='tight')

plt.show()
