import sys
sys.exit(0)

# load default parameter file
cfg.initialize()

# create working directory
wdir = utils.gettempdir('VAS_HEF_wdir')
if not os.path.exists(wdir):
    os.makedirs(wdir)
shutil.rmtree(wdir)
os.makedirs(wdir)
# set path to working directory
cfg.PATHS['working_dir'] = wdir
# set RGI verion and region
cfg.PARAMS['rgi_version'] = rgi_version
# define how many grid points to use around the glacier,
# if you expect the glacier to grow large use a larger border
cfg.PARAMS['border'] = 100
# we use HistAlp climate data
cfg.PARAMS['baseline_climate'] = 'HISTALP'
# set the mb hyper parameters accordingly
cfg.PARAMS['prcp_scaling_factor'] = 1.75
cfg.PARAMS['temp_melt'] = -1.75





# read RGI entry for Hintereisferner as DataFrame
# containing the outline area as shapefile
entity = utils.get_rgi_glacier_entities([rgi_ids]).iloc[0]
# get and set path to intersect shapefile
intersects_db = utils.get_rgi_intersects_region_file(region=rgi_region)
cfg.set_intersects_db(intersects_db)
# initialize the GlacierDirectory
gdir = oggm.GlacierDirectory(entity, reset=True)





# define the local grid and glacier mask
gis.define_glacier_region(gdir, entity=entity)
gis.glacier_masks(gdir)





# run center line preprocessing tasks
centerlines.compute_centerlines(gdir)
centerlines.initialize_flowlines(gdir)
centerlines.compute_downstream_line(gdir)
centerlines.compute_downstream_bedshape(gdir)
centerlines.catchment_area(gdir)
centerlines.catchment_intersections(gdir)
centerlines.catchment_width_geom(gdir)
centerlines.catchment_width_correction(gdir)





# compute the mass balance parameters for the OGGM model
climate.process_histalp_data(gdir)
climate.local_t_star(gdir)
climate.mu_star_calibration(gdir)





# run inversion tasks
inversion.prepare_for_inversion(gdir)
inversion.mass_conservation_inversion(gdir)
inversion.filter_inversion_output(gdir)





# final task
flowline.init_present_time_glacier(gdir)


# In[ ]:


# define number of years to run
nyears = 300
# get equilibrium year t* for HEF
t_star_oggm = gdir.read_json('local_mustar')['t_star']
# the bias is defined to be zero during the calibration process,
# which is why we don't use it here to reproduce the results
cfg.PARAMS['use_bias_for_run'] = False
# run RandomMassBalance model centered around t*, once without tempertaure bias
# and once with positive and negitave temperature bias of 0.5 °C each.
flowline.run_random_climate(gdir, nyears=nyears, y0=t_star_oggm, seed=12, output_filesuffix='_oggm')
flowline.run_random_climate(gdir, nyears=nyears, y0=t_star_oggm, temperature_bias=+0.5,
                             seed=12, output_filesuffix='_oggm_bias_p')
flowline.run_random_climate(gdir, nyears=nyears, y0=t_star_oggm, temperature_bias=-0.5,
                             seed=12, output_filesuffix='_oggm_bias_n');





# define number of years to run
nyears = 300
# get equilibrium year t* for HEF
t_star_oggm = gdir.read_json('local_mustar')['t_star']
# the bias is defined to be zero during the calibration process,
# which is why we don't use it here to reproduce the results
cfg.PARAMS['use_bias_for_run'] = True
# run RandomMassBalance model centered around t*, once without tempertaure bias
# and once with positive and negitave temperature bias of 0.5 °C each.
flowline.run_random_climate(gdir, nyears=nyears, y0=t_star_oggm, seed=12, output_filesuffix='_oggm_mbbias')
flowline.run_random_climate(gdir, nyears=nyears, y0=t_star_oggm, temperature_bias=+0.5,
                             seed=12, output_filesuffix='_oggm_mbbias_bias_p')
flowline.run_random_climate(gdir, nyears=nyears, y0=t_star_oggm, temperature_bias=-0.5,
                             seed=12, output_filesuffix='_oggm_mbbias_bias_n');


# In[ ]:


ds_oggm = normalize_with_start(utils.compile_run_output([gdir], filesuffix='_oggm'))
ds_oggm_p = normalize_with_start(utils.compile_run_output([gdir], filesuffix='_oggm_bias_p'))
ds_oggm_m = normalize_with_start(utils.compile_run_output([gdir], filesuffix='_oggm_bias_n'))


# In[ ]:


ds_oggm_mbbias = normalize_with_start(utils.compile_run_output([gdir], filesuffix='_oggm_mbbias'))
ds_oggm_mbbias_p = normalize_with_start(utils.compile_run_output([gdir], filesuffix='_oggm_mbbias_bias_p'))
ds_oggm_mbbias_n = normalize_with_start(utils.compile_run_output([gdir], filesuffix='_oggm_mbbias_bias_n'))


# In[ ]:


# create figure and axes
fig, [ax0, ax1] = plt.subplots(2, 1, figsize=[6,8])
# plot the evolution of glacier volume
ax0.plot(ds.volume, label='[{} - {}]'.format(t_star-15, t_star+15), c='#2e3131')
ax0.plot(ds_p.volume, label='+0.5 °C', c='#f22613')
ax0.plot(ds_m.volume, label='-0.5 °C', c='#1f3a93')
ax0.axhline(ds.volume[0], c='k', ls=':', lw=0.8, label='initial volume')
ax0.set_ylabel('Relative volume')
ax0.set_title('Volume/Area scaling model')
ax0.legend()
ax1.plot(ds_oggm.volume, label='[{} - {}]'.format(t_star_oggm-15, t_star_oggm+15), c='#6c7a89')
ax1.plot(ds_oggm_p.volume, label='+0.5 °C', c='#f89406')
ax1.plot(ds_oggm_m.volume, label='-0.5 °C', c='#2c82c9')
ax1.axhline(ds_oggm.volume[0], c='k', ls=':', lw=0.8, label='initial volume')
ax1.set_ylabel('Volume [km$^3$]')
ax1.set_title('Relative volume')
ax1.legend()

fig.savefig('vas_oggm_random.pdf', bbox_inches='tight')


# In[ ]:


# create figure and axes
fig, ax0 = plt.subplots(1, 1, figsize=[8,4])

# plot normalized vas volume
ax0.plot(ds.volume, label='[{} - {}]'.format(t_star-15, t_star+15), c='#1f3a93', lw=2)
ax0.plot(ds_p.volume, label='+0.5 °C', c='#4d13d1', lw=2)
ax0.plot(ds_m.volume, label='-0.5 °C', c='#19b5fe', lw=2)

# add aux line at 1
ax0.axhline(ds.volume[0], c='k', ls=':', lw=0.8)

# plot normalized oggm volume
ax0.plot(ds_oggm.volume, label='[{} - {}]'.format(t_star_oggm-15, t_star_oggm+15), c='#d35400', ls='-.', lw=1.5)
ax0.plot(ds_oggm_p.volume, label='+0.5 °C', c='#cf000f', ls='-.', lw=1.5)
ax0.plot(ds_oggm_m.volume, label='-0.5 °C', c='#fcd670', ls='-.', lw=1.5)

# add legend(s)
handels, labels = ax0.get_legend_handles_labels()
l_vas = ax0.legend(handels[:3], labels[:3], title='VAS model', bbox_to_anchor=(1, 1), loc=2)
l_vas.get_title().set_fontweight('bold')
l_oggm = ax0.legend(handels[3:], labels[3:], title='OGGM', bbox_to_anchor=(1, 0), loc=3)
l_oggm.get_title().set_fontweight('bold')
ax0.add_artist(l_vas)

# title, labels, ...
ax0.set_title('Evolution of HEF volume under random (equilibrium) climate')
ax0.set_xlabel('Years of evolution')
ax0.set_ylabel('Relative glacier volume')

# store to file
fig.savefig('vas_oggm_random_oneplot.jpg', bbox_inches='tight')


# In[ ]:


# create figure and axes
fig, ax0 = plt.subplots(1, 1, figsize=[8,4])

# plot normalized vas length
ax0.plot(ds.length, label='[{} - {}]'.format(t_star-15, t_star+15), c='#1f3a93', lw=2)
ax0.plot(ds_p.length, label='+0.5 °C', c='#4d13d1', lw=2)
ax0.plot(ds_m.length, label='-0.5 °C', c='#19b5fe', lw=2)

# add aux line at 1
ax0.axhline(ds.length[0], c='k', ls=':', lw=0.8)

# plot normalized oggm length
ax0.plot(ds_oggm.length, label='[{} - {}]'.format(t_star_oggm-15, t_star_oggm+15), c='#d35400', ls='-.', lw=1.5)
ax0.plot(ds_oggm_p.length, label='+0.5 °C', c='#cf000f', ls='-.', lw=1.5)
ax0.plot(ds_oggm_m.length, label='-0.5 °C', c='#fcd670', ls='-.', lw=1.5)

# add legend(s)
handels, labels = ax0.get_legend_handles_labels()
l_vas = ax0.legend(handels[:3], labels[:3], title='VAS model', bbox_to_anchor=(1, 1), loc=2)
l_vas.get_title().set_fontweight('bold')
l_oggm = ax0.legend(handels[3:], labels[3:], title='OGGM', bbox_to_anchor=(1, 0), loc=3)
l_oggm.get_title().set_fontweight('bold')
ax0.add_artist(l_vas)

# title, labels, ...
ax0.set_title('Evolution of HEF length under random (equilibrium) climate')
ax0.set_xlabel('Years of evolution')
ax0.set_ylabel('Relative glacier length')

# store to file
fig.savefig('vas_oggm_random_oneplot.jpg', bbox_inches='tight')


# ### Power spectrum

# In[ ]:


# load module for the spectral analysis
from scipy import signal


# In[ ]:


# the bias is defined to be zero during the calibration process,
# which is why we don't use it here to reproduce the results
cfg.PARAMS['use_bias_for_run'] = False
# run RandomMassBalance model centered around t*, once without tempertaure bias
# and once with positive and negitave temperature bias of 0.5 °C each.
flowline.run_random_climate(gdir, nyears=3000, y0=t_star_oggm, seed=12, output_filesuffix='_oggm')





# visualize the results
plt.figure(figsize=(5, 4))

# compute the power of the signel per frequency band
sig = ds_lo.length.values.flatten()
freqs, psd = signal.welch(sig)
plt.loglog(freqs, psd, label='vas', c='#FFAE03')

# compute the power of the signel per frequency band
sig = ds_oggm.length.values.flatten()
freqs, psd = signal.welch(sig)
plt.loglog(freqs, psd, label='flowline', c='#004468')

plt.title('PSD of HEF relative length change')
plt.xlabel('Frequency (per year)')
plt.ylabel('Power')
plt.tight_layout()
plt.grid(which='both')
plt.legend()


# ## Finding a historic start area
# TODO...




res = vascaling.find_start_area(gdir, year_start=t_star)
res


# In[ ]:


workflow.execute_entity_task(vascaling.run_random_climate, gdirs, nyears=nyears, y0=t_star, seed=12, init_area_m2=res.x, output_filesuffix='_init')


# In[ ]:


ds_init = utils.compile_run_output([gdir], filesuffix='_init')


# In[ ]:


ax = plt.subplot(1,1,1)
ax.plot(ds.area, label='')
ax.plot(ds_init.area)


# In[ ]:


ds.volume.plot()

