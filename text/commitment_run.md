# VAS commitment run

First I'll deconstruct the commitment run with the Rofental glaciers, given in the [OGGM documentation](http://docs.oggm.org/en/latest/run.html). This includes the `run_random_climate()` task from the `oggm.climate` module.

Next I'll try to figure how to apply it to the VAS model.

## Preprocess a subset of an RGI region

This example shows how to run the first steps of the OGGM preprocessing chain for a subset of the Alps - the Rofental catchment in the Austrian Alps (see [docs.oggm.org](http://docs.oggm.org/en/latest/run_examples/run_rgi_region.html)).

```python
# Python imports
import os
# Libs
import geopandas as gpd
import shapely.geometry as shpg
# Locals
import oggm.cfg as cfg
from oggm import utils, workflow
```

Import external libraries and `oggm` modules.

```python
# Initialize OGGM and set up the default run parameters
cfg.initialize()
rgi_version = '61'
rgi_region = '11'  # Alps

# Local working directory (where OGGM will write its output)
WORKING_DIR = utils.gettempdir('OGGM_Rofental')
utils.mkdir(WORKING_DIR, reset=True)
cfg.PATHS['working_dir'] = WORKING_DIR
```

Load default parameter file, choose the RGI version and region (where 11 corresponds to the Alps). Next step is to create and set the working directory.

```python
# We use intersects
path = utils.get_rgi_intersects_region_file(rgi_region, version=rgi_version)
cfg.set_intersects_db(path)

# RGI file
path = utils.get_rgi_region_file(rgi_region, version=rgi_version)
rgidf = gpd.read_file(path)

# Get the Rofental Basin file
path = utils.get_demo_file('rofental_hydrosheds.shp')
basin = gpd.read_file(path)
```

Get the needed files, which includes the RGI intersects, the RGI entries (including the shapefiles) and the shapefile defining the Rofental.

```python
# Take all glaciers in the Rofental Basin
in_bas = [basin.geometry.contains(shpg.Point(x, y))[0] for
          (x, y) in zip(rgidf.CenLon, rgidf.CenLat)]
rgidf = rgidf.loc[in_bas]
# Store them for later
rgidf.to_file(os.path.join(WORKING_DIR, 'rgi_rofental.shp'))
```

Get all glaciers which center location lies within the Rofental basin.

```python
# Sort for more efficient parallel computing
rgidf = rgidf.sort_values('Area', ascending=False)
cfg.PARAMS['use_multiprocessing'] = True

# Go - initialize glacier directories
gdirs = workflow.init_glacier_regions(rgidf)
```

Since the following tasks are parallelisable, the glaciers are sorted by area to allow for the most efficient computation. As always, the first task  `init_glacier_regions()` initialises the `oggm.GlacierDirections`.

```python
# Tasks shortcuts - see the next examples for more details
workflow.gis_prepro_tasks(gdirs)
workflow.climate_tasks(gdirs)
workflow.inversion_tasks(gdirs)
```

The `workflow` module facilitates the workflow, by combining certain *entity* tasks into thematic packets. The `gis_prepro_tasks` include the tasks computing the glacier grid, computing the centerlines and flowlines, computing the downstream parameters as well as computing the catchment parameters. @TODO: other tasks...

```python
# Compile output
print('Compiling output')
utils.compile_glacier_statistics(gdirs)
utils.write_centerlines_to_shape(gdirs)
```

Lastly, the results are compiled...

## Run random climate

```python
def run_random_climate(gdir, nyears=1000, y0=None, halfsize=15,
                       bias=None, seed=None, temperature_bias=None,
                       store_monthly_step=False,
                       climate_filename='climate_monthly',
                       climate_input_filesuffix='',
                       output_filesuffix='', init_model_fls=None,
                       zero_initial_glacier=False,
                       unique_samples=False,
                       **kwargs):
    """Runs the random mass-balance model for a given number of years.

    This will initialize a
    :py:class:`oggm.core.massbalance.MultipleFlowlineMassBalance`,
    and run a :py:func:`oggm.core.flowline.robust_model_run`.

    Parameters
    ----------
    gdir : :py:class:`oggm.GlacierDirectory`
        the glacier directory to process
    nyears : int
        length of the simulation
    y0 : int, optional
        central year of the random climate period. The default is to be
        centred on t*.
    halfsize : int, optional
        the half-size of the time window (window size = 2 * halfsize + 1)
    bias : float
        bias of the mb model. Default is to use the calibrated one, which
        is often a better idea. For t* experiments it can be useful to set it
        to zero
    seed : int
        seed for the random generator. If you ignore this, the runs will be
        different each time. Setting it to a fixed seed accross glaciers can
        be usefull if you want to have the same climate years for all of them
    temperature_bias : float
        add a bias to the temperature timeseries
    store_monthly_step : bool
        whether to store the diagnostic data at a monthly time step or not
        (default is yearly)
    climate_filename : str
        name of the climate file, e.g. 'climate_monthly' (default) or
        'gcm_data'
    climate_input_filesuffix: str
        filesuffix for the input climate file
    output_filesuffix : str
        this add a suffix to the output file (useful to avoid overwriting
        previous experiments)
    init_model_fls : []
        list of flowlines to use to initialise the model (the default is the
        present_time_glacier file from the glacier directory)
    zero_initial_glacier : bool
        if true, the ice thickness is set to zero before the simulation
    unique_samples: bool
        if true, chosen random mass-balance years will only be available once
        per random climate period-length
        if false, every model year will be chosen from the random climate
        period with the same probability
    kwargs : dict
        kwargs to pass to the FluxBasedModel instance
    """

    mb = MultipleFlowlineMassBalance(gdir, mb_model_class=RandomMassBalance,
                                     y0=y0, halfsize=halfsize,
                                     bias=bias, seed=seed,
                                     filename=climate_filename,
                                     input_filesuffix=climate_input_filesuffix,
                                     unique_samples=unique_samples)

    if temperature_bias is not None:
        mb.temp_bias = temperature_bias

    return robust_model_run(gdir, output_filesuffix=output_filesuffix,
                            mb_model=mb, ys=0, ye=nyears,
                            store_monthly_step=store_monthly_step,
                            init_model_fls=init_model_fls,
                            zero_initial_glacier=zero_initial_glacier,
                            **kwargs)
```