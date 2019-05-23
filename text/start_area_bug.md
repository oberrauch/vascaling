# Finding the best start area 

It's not a bug. It's a feature. Well, kinda...

## What happened so far...

I tried to implement the **V**olume-**A**rea-**S**caling (VAS) model from [Marzeion, et. al. (2012)][] into the OGGM framework. The original model runs with historic climate data (CRU), which calls for a initialisation with a historic glacier area. See following paragraph from the paper.

> Since the relaxation time scale $\tau_A$ introduces memory of past changes into the model, it is not possible to integrate the model backwards in time to determine the evolution of the glacier before the year of surface area measurement. For this reason, the glacier’s surface area $A_{start}$ at the beginning of the model integration (i.e., 1901 for the forcing with observed climate variability and change, and 1850 for most cases of modeled climate variability and change) is estimated by iteratively seeking that surface area in the starting year of the integration that will result in the measured surface area in the year of the measurement. The iteration is deemed successful when the modeled surface area is within 0.1% of the measured surface area during the year of the measurement; the iterative process is broken off after 100 iterations if unsuccessful (see Sect. 6.2.2 how these glaciers are treated).
>
> [Marzeion, et. al. (2012)][]

So, after the mass balance model and the "dynamic" scaling model worked, I tried to implement the process of seeking the start area for a certain starting year in the past. However, I use an optimisation (minimisation) function rather than an iterative process. I'll explain my code and the reasoning behind it in the following section.

### My code

```python
def find_start_area(gdir, year_start=1851):
    """ This task find the start area for the given glacier, which results in
    the best results after the model integration (i.e., modeled glacier surface
    closest to measured RGI surface in 2003).

    All necessary prepro task (gis, centerline, climate) must be executed
    beforehand, as well as the local_t_star() task.

    :param gdir: (oggm.GlacierDirectory)
    :param year_start: (int, optional) year at the beginning of the model
        integration, default = 1851 (best if working with HISTALP)
    :return: (scipy.optimize.OptimizeResult)
    """
```

The tasks `find_start_area()` need a `oggm.GalcierDirectory gdir` as input. Additionally it is possible to define the year `year_start` at the beginning of the model integration. The default year is 1851, since I'm mostly working with HistAlp data. All needed preprocessing tasks (i.e., GIS, centerline, climate tasks) must be executed beforehand, in order to allow the initialization of a mass balance model from the `gdir`.

```python
    # instance the mass balance models
    mbmod = VAScalingMassBalance(gdir)

    # get reference area and year from RGI
    a_rgi = gdir.rgi_area_m2
    rgi_df = utils.get_rgi_glacier_entities([gdir.rgi_id])
    y_rgi = int(rgi_df.BgnDate.values[0][:4])
    # get min and max glacier surface elevation
    h_min, h_max = get_min_max_elevation(gdir)

    # set up the glacier model with the reference values (from RGI)
    model_ref = VAScalingModel(year_0=y_rgi, area_m2_0=a_rgi,
                               min_hgt=h_min, max_hgt=h_max,
                               mb_model=mbmod)
```



```python
    def _to_minimize(area_m2_start, ref, year_start=year_start):
        """ Initialize VAS glacier model as copy of the reference model (ref)
        and adjust the model to the given starting area (area_m2_start) and
        starting year (1851). Let the model evolve to the same year as the
        reference model. Compute and return the relative absolute area error.

        :param area_m2_start: (float)
        :param ref: (oggm.VAScalingModel)
        :return: (float) relative absolute area error
        """
        # define model
        model_tmp = VAScalingModel(year_0=ref.year_0,
                                   area_m2_0=ref.area_m2_0,
                                   min_hgt=ref.min_hgt_0,
                                   max_hgt=ref.max_hgt,
                                   mb_model=ref.mb_model)
        # scale to desired starting size
        model_tmp.create_start_glacier(area_m2_start, year_start=year_start)
        # run and compare, return relative error
        return np.abs(model_tmp.run_and_compare(ref))

    # define bounds - between 100m2 and two times the reference size
    area_m2_bounds = [100, 2 * model_ref.area_m2_0]
    # run minimization
    minimization_res = minimize_scalar(_to_minimize, args=(model_ref),
                                       bounds=area_m2_bounds,
                                       method='bounded')

    return minimization_res
```



## References

[Marzeion, et. al. (2012)]: http://www.the-cryosphere.net/6/1295/2012/	"Past and future sea-level change from the surface mass balance of glaciers"