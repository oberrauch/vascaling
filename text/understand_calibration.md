# OGGM $\mu^*$ calibration

Hereafter I try to understand the code used for the calibration of the temperature sensitivity $\mu^*$ in OGGM. The calibration process is based on two tasks performed on the GlacierDirectory, `local_t_start(gdir)` and `mu_star_calibration(gdir)`. 

The  `local_t_start(gdir)`  task computes the local $t^*$ (simplified the year around which the glacier was in equilibrium) and the associated glacier wide $\mu^*$ for the given glacier. The function adds the relevant mass balance parameters to the `climate_info.pkl` file, the parameters being the default temperature ggradiet $\gamma_\text{temp}$, the temperature below which all precipitation falls in form of snow $T^\text{prec solid}$, the temperature above which all precipitation falls in form of rain $T^\text{prec liquid}$, the temperature above which ice melt occurs $T^\text{melt}$, as well as the precipitation scaling factor $a$. This allows other tasks to access this information and at the same time inhibits other accidental changes on the calibration parameters. The actual output of the function is written into the `local_mustar.json` file, containing the glacier's RGI ID, the year (or index) $t^*$, the bias $\beta^*$ and the glacier wide $\mu^*$.

The `mu_star_calibration(gdir)` calibrates the temperature sensitivity $\mu^*$ and computes the apparent mass balance for each flow line. Apparent mass balance refers to the average specific mass balance for the climatological period around $t^*$, whereby the specific mass balance is the mass balance averaged over the glaciers surface ([IPPC, WG 1](https://www.ipcc.ch/ipccreports/tar/wg1/413.htm)). This is done via a recursive process, which takes care of eventually arising negative mass fluxes on the last grid point(s). The flow lines are stored in the `inversion_flowlines.pkl` file, each with its own $\mu^*$. Additional parameters are written to the `local_mustat.json` file, containing the $\mu^*$  of each flow line, a flow line averaged $\overline{\mu^*}$ as well as a flag whether or not the $\mu^*$ is the same for all flow lines.

## Theoretical background

Up next, a short summary of the main calibration ideas developed by [Marzeion et. al., 2010](http://www.marzeion.info/sites/default/files/marzeion_etal_12a.pdf). 

**Mass balance model**

Simplified, a glaciers (annual specific surface) mass balance $B$ is the difference between accumulation (i.e. mass gain by snowfall, avalanches, drift, ...) and ablation (i.e. mass loss via ice melt, sublimation, ...) over the course of a year. The used mass balance model relies solely on the area mean monthly solid precipitation onto the glacier surface $P_i^{\text{solid}}$ and the monthly mean air temperature at the location and elevation of the terminus of the glacier $T_i^{\text{terminus}}$ as input. 
$$
B = \left[\sum_{i=1}^{12}\left[
			P_{i}^{\text{solid}} 
			- \mu^* \cdot
					\max\left(T_{i}^{\text{terminus}}
					- T_{\text{melt}},\ 0\right)
	\right]\right] - \beta^*
$$
Hereby, $\mu^*$ is the temperature sensitivity of the glacier and $T_{\text{melt}}$ the mean air temperature above which ice melt occurs.

**Temperature sensitivity for glaciers with mass balance measurements**

The first step is to estimate the temperature sensitivity $\mu^*$ for all glaciers with available (and reliable) mass balance measurements, which is a global total of 255 glaciers used by [Marzeion et. al., 2010](http://www.marzeion.info/sites/default/files/marzeion_etal_12a.pdf). For all this glaciers a *yearly* temperature sensitivity $\mu(t)$ is calculated, whereby the year $t$ indicates the centre of a 31 year climatological period, assuming that the glacier was in equilibrium (i.e. mass balance $B = 0$) during that period. In order words, solving the following equation for $\mu(t)$.
$$
B = \sum_{i=1}^{12}\left[
			P(t)_{i,\text{ clim}}^{\text{solid}} 
			- \mu(t) \cdot
					\max\left(T(t)_{i,\text{ clim}}^{\text{terminus}}
					- T_{\text{melt}},\ 0\right)
	\right] = 0
$$
Hereby, $P(t)_{i,\text{ clim}}^{\text{solid}}$ and $T(t)_{i,\text{ clim}}^{\text{terminus}}$ are the climatological values of the area mean monthly solid precipitation onto the glacier surface $P_i^{\text{solid}}$ and the monthly mean air temperature at the location and elevation of the terminus of the glacier $T_i^{\text{terminus}}$, respectively. This yields a number of values of $\mu(t)$ for each glacier,  whereby the year $t$ can be understood as an index for different climate conditions. For each glacier, each of its computed temperature sensitivities $\mu(t)$ is than used to calculate the mass balance for the years of mass balance, using above defined mass balance model equation. The goal is to identify a year $t^*$ (or much rather the climate conditions around said year) for which the absolute difference $\beta^*$ between the mean of the modelled mass balances during the years of mass balance measurements $\overline{B(t)_{\text{modelled}}}$ and the mean of the observed mass balances $\overline{B_{\text{measured}}}$ is minimal, i.e.
$$
\left|\overline{B(t)_{\text{modelled}}} - \overline{B_{\text{measured}}}\right| = \beta^*.
$$
In order words, for each glacier with mass balance measurements a temperature sensitivity $\mu^* = \mu(t^*)$ is determined that produces the smallest possible bias $\beta^* = \beta(t^*)$. Thereby, the year $t^*$ is the index to the climate conditions which produces the best model results. This does not necessarily mean, that the glacier was in equilibrium around $t^*$.

**Temperature sensitivity for glaciers without mass balance measurements**

The "revolutionary" idea of [Marzeion et. al., 2010](http://www.marzeion.info/sites/default/files/marzeion_etal_12a.pdf) was to interpolate $t^*$ to glaciers without mass balance measurements, much rather than interpolating the temperature sensitivity $\mu^*$. The interpolation is done by using the ten nearest glaciers with mass balance measurements, weighting inversely with distance. Given a $t^*$, the temperature sensitivity $\mu^*$ is determined by requiring that
$$
B = \left[\sum_{i=1}^{12}\left[
			P(t^*)_{i,\text{ clim}}^{\text{solid}} 
			- \mu^* \cdot
					\max\left(T(t^*)_{i,\text{ clim}}^{\text{terminus}}
					- T_{\text{melt}},\ 0\right)
	\right]\right] = 0.
$$
In words: The temperature sensitivity $\mu^*$ is chosen, so that the calculated mass balance equals zero, using the monthly climatological values around $t^*$ of solid precipitation and temperature as input,  $P(t^*)_{i,\text{ clim}}^{\text{solid}}$ and $T(t^*)_{i,\text{ clim}}^{\text{terminus}}$ respectively. In even simpler words: solve the above equation for $\mu^*$.

**Bias correction**

@TODO: where is the bias used in OGGM

## `local_t_star(gdir)`

`local_t_star(gdir)` is the first task to be executed. It computes the local $t^*$ (simplified the year around which the glacier was in equilibrium) and the associated glacier wide $\mu^*$ for the given glacier.

```python
@entity_task(log, writes=['local_mustar'])
def local_t_star(gdir, *, ref_df=None, tstar=None, bias=None):
    """Compute the local t* and associated glacier-wide mu*.

    If ``tstar`` and ``bias`` are not provided, they will be interpolated from
    the reference t* list.

    Note: the glacier wide mu* is here just for indication. It might be
    different from the flowlines' mu* in some cases.

    Parameters
    ----------
    gdir : oggm.GlacierDirectory
    ref_df : pd.Dataframe, optional
        replace the default calibration list with your own.
    tstar: int, optional
        the year where the glacier should be equilibrium
    bias: float, optional
        the associated reference bias
    """
```

The task operates on a GlacierDirectory `gdir` , specifying the glacier. It is possible to supply costume values for  $t^*$ ( `tstar`) and the bias $\beta^*$( `bias`), otherwise they will be interpolated from the reference list. Thereby it is also possible to provide a costume calibration list `ref_df`.

**Define relevant mass balance calibration parameters**

```python
    # Relevant mb params
    params = ['temp_default_gradient', 'temp_all_solid', 'temp_all_liq',
              'temp_melt', 'prcp_scaling_factor']
```

The mass balance parameters are later written in the `climate_info.pkl` file, which stores necessary climate parameters for other task. This approach makes sure that no other task messes with those parameters:

- `temp_default_gradient` $\gamma_\text{temp}$ is the default temperature gradient in Kelvin per m.
- `temp_all_solid` $T^\text{prec solid}$ is the temperature below which all precipitation falls as snow
- `temp_all_liq` $T^\text{prec liquid}$ is the temperature above which all precipitation falls as rain
- `temp_melt` $T^\text{melt}$ is the temperature above which the glacier ice starts to melt
- `prcp_scaling_factor` $a$ multiplies the given precipitation amount in the climate file

**Interpolation of $t^*$**

```python
    if tstar is None or bias is None:
        # Do our own interpolation
        if ref_df is None:
            if not cfg.PARAMS['run_mb_calibration']:
                # Make some checks and use the default one
                climate_info = gdir.read_pickle('climate_info')
                source = climate_info['baseline_climate_source']
                ok_source = ['CRU TS4.01', 'CRU TS3.23', 'HISTALP']
                if not np.any(s in source.upper() for s in ok_source):
                    msg = ('If you are using a custom climate file you should '
                           'run your own MB calibration.')
                    raise MassBalanceCalibrationError(msg)
                v = gdir.rgi_version[0]  # major version relevant

                # Check that the params are fine
                str_s = 'cru4' if 'CRU' in source else 'histalp'
                vn = 'ref_tstars_rgi{}_{}_calib_params'.format(v, str_s)
                for k in params:
                    if cfg.PARAMS[k] != cfg.PARAMS[vn][k]:
                        raise ValueError('The reference t* you are trying '
                                         'to use was calibrated with '
                                         'difference MB parameters. You '
                                         'might have to run the calibration '
                                         'manually.')
                ref_df = cfg.PARAMS['ref_tstars_rgi{}_{}'.format(v, str_s)]
            else:
                # Use the the local calibration
                fp = os.path.join(cfg.PATHS['working_dir'], 'ref_tstars.csv')
                ref_df = pd.read_csv(fp)
```
If $t^*$ and $\beta^*$ are not provided, they will be calibrated from the reference list including all glaciers with mass balance measurements. It is possible to supply a costume reference list via `ref_df`. If a costume climate file (i.e. not *CRU* or *HistAlp*) is used the“ mass balance calibration should be run and the local calibration file should be used. Otherwise the default reference list is used. Hereby the parameters are compared to the default parameters for the given combination of RGI version and climate file.

```python
        # Compute the distance to each glacier
        distances = utils.haversine(gdir.cenlon, gdir.cenlat,
                                    ref_df.lon, ref_df.lat)

        # Take the 10 closest
        aso = np.argsort(distances)[0:9]
        amin = ref_df.iloc[aso]
        distances = distances[aso]**2

        # If really close no need to divide, else weighted average
        if distances.iloc[0] <= 0.1:
            tstar = amin.tstar.iloc[0]
            bias = amin.bias.iloc[0]
        else:
            tstar = int(np.average(amin.tstar, weights=1./distances))
            bias = np.average(amin.bias, weights=1./distances)
```
After a reference list `ref_df` is chosen, the actual interpolation can start. First the distance between the given glacier and all glaciers of the reference list is computed. The function `haversine(lon1, lat1, lon2, lat2)` from the `utils` module computes the great circle distance between two (or more) points on Earth. The $t^*$ value of ten closest glaciers with mass balance measurements is interpolated, weighting inversely with distance. Given a very close glacier (i.e. with a distance of less than or equal to 0.1 meters) its value is used and no interpolation is needed.

**Mass balance parameters**

```python
    # Add the climate related params to the GlacierDir to make sure
    # other tools cannot fool around without re-calibration
    out = gdir.read_pickle('climate_info')
    out['mb_calib_params'] = {k: cfg.PARAMS[k] for k in params}
    gdir.write_pickle(out, 'climate_info')
```
Here the mass balance calibration parameters are written into the `climate_info.pkl` file, as mentioned above.

**Computing $\mu^*$**

```python
    # We compute the overall mu* here but this is mostly for testing
    # Climate period
    mu_hp = int(cfg.PARAMS['mu_star_halfperiod'])
    yr = [tstar - mu_hp, tstar + mu_hp]

    # Do we have a calving glacier?
    cmb = calving_mb(gdir)
```
The climate period is defined as 31 years centred around $t^*$ (if `mu_star_halfperiod` is not altered in the `PARAMS.cfg` file). Calving has to be taken into account since it drastically alters the mass balance.

```python
    # Get the corresponding mu
    years, temp_yr, prcp_yr = mb_yearly_climate_on_glacier(gdir, year_range=yr)
    assert len(years) == (2 * mu_hp + 1)
```
In order to compute $\mu^*$ the meteorological parameters relevant for the mass balance (i.e. positive temperature and solid precipitation) during the climate period need to be know. This information is gathered using the `mb_yearly_climate_on_glacier(gdir)` function, which is explained in more detail below.

```python
    # mustar is taking calving into account  (units of specific MB)
    mustar = (np.mean(prcp_yr) - cmb) / np.mean(temp_yr)
```

This is the actual calculation, based on the equation described in the theoretical background.

```python
	if not np.isfinite(mustar):
            raise MassBalanceCalibrationError('{} has a non finite '
                                              'mu'.format(gdir.rgi_id))
	# Clip the mu
    if not (cfg.PARAMS['min_mu_star'] < mustar < cfg.PARAMS['max_mu_star']):
    	raise MassBalanceCalibrationError('mu* out of specified bounds.')
```

The calculation can yield in a non physical result, for which is tested here.

```python
    # Scalars in a small dict for later
    df = dict()
    df['rgi_id'] = gdir.rgi_id
    df['t_star'] = int(tstar)
    df['bias'] = bias
    df['mu_star_glacierwide'] = mustar
    gdir.write_json(df, 'local_mustar')
```
The result of this calibration/interpolation/computation process is saved in the `local_mustar.json` file, making the following parameters accessible to other tasks:

- `rgi_id` is the RGI ID in order to identify the glacier
- `t_start` is most obviously $t^*$
- `bias` is most obviously the bias $\beta^*$
- `mu_star_glacierwide` is the glacier wide temperature sensitivity $\mu^*$



### `mb_yearly_climate_on_glacier(gdir)`

This function is called by the `locat_t_star(gdir)` task and computes the so called glacier wide *mass balance climate*, i.e. the meteorological parameters relevant for mass gain and loss. It returns the temperature energies (causing melt) and solid precipitation for each year and each elevation band. To do so it calls the `mb_yearly_climate_on_height(gdir, heights)` function, which itself uses the `mb_climate_on_height(gdir, heights)` function. The latter actually reads the glaciers climate file.

```python
def mb_yearly_climate_on_glacier(gdir, *, year_range=None):
    """Yearly mass-balance climate at all glacier heights,
    multiplied with the flowlines widths. (all in pix coords.)

    The precipitation time series are not corrected! @ASK

    Parameters
    ----------
    gdir : GlacierDirectory
        the glacier directory
    year_range : [int, int], optional
        Provide a [y0, y1] year range to get the data for specific
        (hydrological) years only.

    Returns
    -------
    (years, tempformelt, prcpsol)::
        - years: array of shape (ny,)
        - tempformelt:  array of shape (len(heights), ny)
        - prcpsol:  array of shape (len(heights), ny) (not corrected!)
    """
```

As most other tasks it operates on a GlacierDirectory `gdir`.  It is possible to specify a range of (hydrological) years `year_range`, default is to read all available climate data.

```python
    flowlines = gdir.read_pickle('inversion_flowlines')

    heights = np.array([])
    widths = np.array([])
    for fl in flowlines:
        heights = np.append(heights, fl.surface_h)
        widths = np.append(widths, fl.widths)
```
Given the inversion flow line(s) all elevation band and corresponding glacier widths are gathered.

```python
    years, temp, prcp = mb_yearly_climate_on_height(gdir, heights,
                                                    year_range=year_range,
                                                    flatten=False)
```
The `mb_yearly_climate_on_height()` function returns positive temperatures and solid precipitation for each elevation band, details below.

```python
    temp = np.average(temp, axis=0, weights=widths)
    prcp = np.average(prcp, axis=0, weights=widths)

    return years, temp, prcp
```

Temperature and precipitation are averaged over the entire glacier, weighted by width of the corresponding elevation band.



### `mb_yearly_climate_on_height(gdir, heights)`

This function uses the underlying `mb_climate_on_height(gdir, heights)` to compute yearly values (sums) of the climate inputs driving the mass balance for a given glacier and for given elevation levels.

```python
def mb_yearly_climate_on_height(gdir, heights, *,
                                year_range=None, flatten=False):
    """Yearly mass-balance climate of the glacier at a specific height

    The precipitation time series are not corrected!

    Parameters
    ----------
    gdir : GlacierDirectory
        the glacier directory
    heights: ndarray
        a 1D array of the heights (in meter) where you want the data
    year_range : [int, int], optional
        Provide a [y0, y1] year range to get the data for specific
        (hydrological) years only.
    flatten : bool
        for some applications (glacier average MB) it's ok to flatten the
        data (average over height) prior to annual summing.

    Returns
    -------
    (years, tempformelt, prcpsol)::
        - years: array of shape (ny,)
        - tempformelt:  array of shape (len(heights), ny) (or ny if flatten
        is set)
        - prcpsol:  array of shape (len(heights), ny) (or ny if flatten
        is set)
    """
```

As all functions above it operates on the GlacierDirectory `gdir`. The elevation levels `heights` must be supplied. It is possible to specify a range of (hydrological) years `year_range`. Default is to read all available climate data. In addition there is a possibility to compute a spatial average `flatten` of the returned parameters. Default is to return a distinct value for each given elevation level.

```python
	# [Get monthly mass balance climate parameters]
	time, temp, prcp = mb_climate_on_height(gdir, heights,
                                            year_range=year_range)

    # [Check if climate data includes all 12 month of all years]
    ny, r = divmod(len(time), 12)
    if r != 0:
        raise ValueError('Climate data should be N full years exclusively')
    # Last year gives the tone of the hydro year
    years = np.arange(time[-1].year-ny+1, time[-1].year+1, 1)
```

The function `mb_climate_on_height(gdir, heights)` returns a two dimensional arrays for temperature above freezing point `temp2dformelt` and solid precipitation `prcpsol` for all elevation levels `heights` and all month in the given range of (hydrological) years `year_range`. See below for details.

After making sure that all years consist of 12 months, the hydrological year is computed from the given timestamps. Hereby the last year is decisive, e.g. the hydrological year 2010 starts on October 1, 2009 and ends on September 30, 2010 ([more on wikipedia.org](https://en.wikipedia.org/wiki/Water_year)).

```python
    if flatten:
        # Spatial average
        temp_yr = np.zeros(len(years))
        prcp_yr = np.zeros(len(years))
        temp = np.mean(temp, axis=0)
        prcp = np.mean(prcp, axis=0)
        for i, y in enumerate(years):
            temp_yr[i] = np.sum(temp[i*12:(i+1)*12])
            prcp_yr[i] = np.sum(prcp[i*12:(i+1)*12])
    else:
        # Annual prcp and temp for each point (no spatial average)
        temp_yr = np.zeros((len(heights), len(years)))
        prcp_yr = np.zeros((len(heights), len(years)))
        for i, y in enumerate(years):
            temp_yr[:, i] = np.sum(temp[:, i*12:(i+1)*12], axis=1)
            prcp_yr[:, i] = np.sum(prcp[:, i*12:(i+1)*12], axis=1)
```
To obtain yearly values the monthly temperature above freezing and the solid precipitation are summed up, which makes sense since we talk about yearly energy and mass input, respectively. The function gives the option to compute a spatial average over the entire glacier for each year, no weighting depending on glacier width at given elevation level. Otherwise (which is default) each elevation level has its own values.

```python
    return years, temp_yr, prcp_yr
```

The function returns the hydrological year `years`, and yearly values of temperature above freezing point `temp_yr` and solid precipitation `prcp_yr`. The shape of those parameters depends on the chosen option on spatial averaging. If spatial averaging is applied all there is only one value of temperature and precipitation per year. Otherwise the parameters come in a matrix shape, whereby the second dimension accounts for the different elevation levels.



### `mb_climate_on_height(gdir, heights)`

This is the *base level* function which actually reads the climate file and computes the climate inputs driving the mass balance for a given glacier and for given elevation levels.

```python
def mb_climate_on_height(gdir, heights, *, time_range=None, year_range=None):
    """Mass-balance climate of the glacier at a specific height

    Reads the glacier's monthly climate data file and computes the
    temperature "energies" (temp above 0) and solid precipitation at the
    required height.

    Parameters
    ----------
    gdir : GlacierDirectory
        the glacier directory
    heights: ndarray
        a 1D array of the heights (in meter) where you want the data
    time_range : [datetime, datetime], optional
        default is to read all data but with this you
        can provide a [t0, t1] bounds (inclusive).
    year_range : [int, int], optional
        Provide a [y0, y1] year range to get the data for specific
        (hydrological) years only. Easier to use than the time bounds above.

    Returns
    -------
    (time, tempformelt, prcpsol)::
        - time: array of shape (nt,)
        - tempformelt:  array of shape (len(heights), nt)
        - prcpsol:  array of shape (len(heights), nt)
    """
```

As all functions above it operates on the GlacierDirectory `gdir`. The elevation levels `heights` must be supplied. It is possible to specify a time range `time_range` by `datetime` objects or a range of (hydrological) years `year_range`. Default is to read all available climate data.

```python
    if year_range is not None:
        sm = cfg.PARAMS['hydro_month_' + gdir.hemisphere]
        em = sm - 1 if (sm > 1) else 12
        t0 = datetime.datetime(year_range[0]-1, sm, 1)
        t1 = datetime.datetime(year_range[1], em, 1)
        return mb_climate_on_height(gdir, heights, time_range=[t0, t1])
```
The possibly supplied year range is converted into `datetime` format accounting for differences in the hydrological cycle in the norther and southern hemisphere.

```python
    # Parameters
    temp_all_solid = cfg.PARAMS['temp_all_solid']
    temp_all_liq = cfg.PARAMS['temp_all_liq']
    temp_melt = cfg.PARAMS['temp_melt']
    prcp_fac = cfg.PARAMS['prcp_scaling_factor']
    default_grad = cfg.PARAMS['temp_default_gradient']
    g_minmax = cfg.PARAMS['temp_local_gradient_bounds']
```

Here the necessary mass balance parameters are gathered from the `PARAMS.cfg` file:

- `temp_all_solid` $T^\text{prec solid}$ is the temperature below which all precipitation falls as snow.
- `temp_all_liq` $T^\text{prec liquid}$ is the temperature above which all precipitation falls as rain.
- `temp_melt` $T^\text{melt}$ is the temperature above which the glacier ice starts to melt.
- `prcp_scaling_factor` $a$ multiplies the given precipitation amount in the climate file.
- `temp_default_gradient` $\gamma_\text{temp}$ is the default temperature gradient in Kelvin per m.
- `g_minmax` holds physical (realistic) boundaries for the temperature gradient.

```python
    # Read file
    igrad = None
    with utils.ncDataset(gdir.get_filepath('climate_monthly'), mode='r') as nc:
        # time
        time = nc.variables['time']
        time = netCDF4.num2date(time[:], time.units)
        if time_range is not None:
            p0 = np.where(time == time_range[0])[0]
            try:
                p0 = p0[0]
            except IndexError:
                raise MassBalanceCalibrationError('time_range[0] not found in '
                                                  'file')
            p1 = np.where(time == time_range[1])[0]
            try:
                p1 = p1[0]
            except IndexError:
                raise MassBalanceCalibrationError('time_range[1] not found in '
                                                  'file')
        else:
            p0 = 0
            p1 = len(time)-1

        time = time[p0:p1+1]

        # Read timeseries
        itemp = nc.variables['temp'][p0:p1+1]
        iprcp = nc.variables['prcp'][p0:p1+1]
        if 'gradient' in nc.variables:
            igrad = nc.variables['gradient'][p0:p1+1]
            # Security for stuff that can happen with local gradients
            igrad = np.where(~np.isfinite(igrad), default_grad, igrad)
            igrad = np.clip(igrad, g_minmax[0], g_minmax[1])
        ref_hgt = nc.ref_hgt

        # Default gradient?
        if igrad is None:
            igrad = itemp * 0 + default_grad

        # Correct precipitation
        iprcp *= prcp_fac
```

Reading the monthly climate file. If a time period is specified, the climate data will be trimmed accordingly given the range is covered. Monthly temperature `itemp` and precipitation `iprcp` are directly adopted, whereby the precipitation amount is corrected using the given scaling factor. The gradient `igrad`, if supplied by the climate file, is corrected for non finite and non physical values, otherwise the default gradient is used. Furthermore the reference height `ref_hgt` is stored.

```python
    # For each height pixel:
    # Compute temp and tempformelt (temperature above melting threshold)
    npix = len(heights)
    grad_temp = np.atleast_2d(igrad).repeat(npix, axis=0)
    grad_temp *= (heights.repeat(len(time)).reshape(grad_temp.shape) - ref_hgt)
    temp2d = np.atleast_2d(itemp).repeat(npix, axis=0) + grad_temp
    temp2dformelt = temp2d - temp_melt
    temp2dformelt = np.clip(temp2dformelt, 0, temp2dformelt.max())
```

The temperature for each elevation level `temp2d` is computed, using the monthly temperature gradient `igrad`, the height levels `heights` and the temperature `itemp` at the reference elevation. In more detail: The gradient is given as an array of length $\text{#(months)}$. In order to compute the temperature for every height pixel it is transformed into a matrix of shape $\text{#(months)}\times\text{#(height levels)}$ .  In a second step, the elevation differences with respect to the reference height are computed (in the same matrix format as the gradient) and multiplied entry wise with the gradient matrix. Finally, adding the reference temperature (brought into the same shape) yields a monthly temperature for each elevation level.

The *temperature energy* `temp2dformelt` results from the difference of actual temperature `temp2d`  and the melting threshold `temp_melt`, whereby only positive values are considered.

```python
    # [For each height pixel:]
    # Compute solid precipitation from total precipitation
    prcpsol = np.atleast_2d(iprcp).repeat(npix, axis=0)
    fac = 1 - (temp2d - temp_all_solid) / (temp_all_liq - temp_all_solid)
    fac = np.clip(fac, 0, 1)
    prcpsol = prcpsol * fac
```

The amount of solid precipitation `prcpsol` at a given altitude is computed as a certain percentage of the overall precipitation amount `iprcp`. Thereby, the scaling factor `fac` depends on the temperature at the respective altitude. The factor is one (i.e., total precipitation amount is solid) below the given lower threshold `temp_all_solid`, linearly decreasing to zero (i.e., total precipitation amount is liquid) at the given upper threshold `temp_all_liquid`.

```python
    return time, temp2dformelt, prcpsol
```

The function returns the date (month/year) `time`, and a two dimensional array of temperature above freezing point `temp2dformelt` and solid precipitation `prcpsol` for all elevation levels and all dates.



## `mu_star_calibration(gdir)`

The `mu_star_calibration(gdir)` function is called after interpolating $t^*$ and computing the corresponding glacier wide $\mu^*$. It calibrates the temperature sensitivity $\mu^*$ for each flow line and computes the associated apparent mass balance.

```python
def mu_star_calibration(gdir):
    """Compute the flowlines' mu* and the associated apparent mass-balance.

    Parameters
    ----------
    gdir : oggm.GlacierDirectory
    """
```

As all other entity task it operates on a GlacierDirectory `gdir`.

```python
    # Interpolated data
    df = gdir.read_json('local_mustar')
    t_star = df['t_star']
    bias = df['bias']

    # For each flowline compute the apparent MB
    fls = gdir.read_pickle('inversion_flowlines')
    # If someone calls the task a second time we need to reset this
    for fl in fls:
        fl.mu_star_is_valid = False
```

The calibration process relies on $t^*$ and $\beta^*$, which are both read from the `local_mustar.json` file. The next step is to load the glacier flow lines from the `inversion_flowlines.pkl` file, whereby the possibly set `mu_star_is_valid` flag is reset to `False`. 

```python
    # Let's go
    _recursive_mu_star_calibration(gdir, fls, t_star)
```

This task does the actual calibration work, and is explained in more detail below. In summary, it calibrates the $\mu^*$ for each flow line. If told to do so via the `correct_for_neg_flux` flag, the function takes care of arising problems with negative fluxes using a recursive process.

**Attention**: there are two different but very similar flags in the `PARAMS.cfg` file, concerning problems with negative fluxes. The `correct_for_neg_flux` flag enables a recursive $\mu^*$ calibration process for each flow line until all computed fluxes are physical (using the above function). Whereas, if the `filter_for_neg_flux` flag is set, all erroneous flow lines are deleted and computed all over (see block below).

```python
    # If the user wants to filter the bad ones we remove them and start all
    # over again until all tributaries are physically consistent with one mu
    do_filter = [fl.flux_needs_correction for fl in fls]
    if cfg.PARAMS['filter_for_neg_flux'] and np.any(do_filter):
        assert not do_filter[-1]  # This should not happen
        # Keep only the good lines
        # TODO: this should use centerline.line_inflows for more efficiency!
        heads = [fl.orig_head for fl in fls if not fl.flux_needs_correction]
        centerlines.compute_centerlines(gdir, heads=heads, reset=True)
        centerlines.initialize_flowlines(gdir, reset=True)
        if gdir.has_file('downstream_line'):
            centerlines.compute_downstream_line(gdir, reset=True)
            centerlines.compute_downstream_bedshape(gdir, reset=True)
        centerlines.catchment_area(gdir, reset=True)
        centerlines.catchment_intersections(gdir, reset=True)
        centerlines.catchment_width_geom(gdir, reset=True)
        centerlines.catchment_width_correction(gdir, reset=True)
        local_t_star(gdir, tstar=t_star, bias=bias, reset=True)
        # Ok, re-call ourselves
        return mu_star_calibration(gdir, reset=True)
```

The flow lines are newly computed starting from the heads of the existing and physically consistent flow lines. This includes the following tasks (which are mostly self explanatory, for more details see `oggm.core.centerlines`):

- `compute_centerlines(gdir)`: compute the center lines following [Kienholz et. al., (2014)](https://www.the-cryosphere.net/8/503/2014/).
- `initialize_flowlines(gdir)`: transforms the geometrical `Centerlines` in the more physical
  ​    *Inversion Flow Lines*.
- `compute_downstream_line(gdir)`: compute the line continuing the glacier.
- `compute_downstream_bedshape(gdir)`: the bed shape obtained by fitting a parabola to the line's normals. Also computes the downstream's altitude.
- `catchment_area(gdir)`: compute the catchment areas of each tributary line.
- `catchment_intersections(gdir)`: computes the intersections between the catchments.
- `catchment_width_geom(gdir)`: compute geometrical catchment widths for each point of the flow lines.
- `catchment_width_correction(gdir)`: corrects for `NaN`s and inconsistencies in the geometrical widths.

Afterwards, the local $t^*$ and the corresponding glacier wide $\mu^*$ are recomputed using `local_t_star(gdir)`. It follows the recursive call of  `mu_star_calibration(gdir)` , which then again calibrates $\mu^*$ for each flow line until all fluxes are positive eventually.

```python
    # Check and write
    rho = cfg.PARAMS['ice_density']
    aflux = fls[-1].flux[-1] * 1e-9 / rho * gdir.grid.dx**2
    # If not marine and a bit far from zero, warning
    cmb = calving_mb(gdir)
    if cmb == 0 and not np.allclose(fls[-1].flux[-1], 0., atol=0.01):
        log.warning('(%s) flux should be zero, but is: '
                    '%.4f km3 ice yr-1', gdir.rgi_id, aflux)
    # If not marine and quite far from zero, error
    if cmb == 0 and not np.allclose(fls[-1].flux[-1], 0., atol=1):
        msg = ('({}) flux should be zero, but is: {:.4f} km3 ice yr-1'
               .format(gdir.rgi_id, aflux))
        raise MassBalanceCalibrationError(msg)
    gdir.write_pickle(fls, 'inversion_flowlines')
```
For this last step of quality control the flux through the last grid point is converted into understandable units of volume ice per year $[\mathrm{km^3\ ice\ yr^{-1}}]$. The flux should be zero, since the glacier is supposably in equilibrium during the climatological period around $t^*$. Small fluxes below $1\ \mathrm{km^3\ ice\ yr^{-1}}$ raise a warning, whereas larger fluxes larger result in an error.

The flow lines with there $\mu^*$ are stored in the `inversion_flowlines.pkl` file.

```python
    # Store diagnostics
    mus = []
    weights = []
    for fl in fls:
        mus.append(fl.mu_star)
        weights.append(np.sum(fl.widths))
    df['mu_star_per_flowline'] = mus
    df['mu_star_flowline_avg'] = np.average(mus, weights=weights)
    all_same = np.allclose(mus, mus[0], atol=1e-3)
    df['mu_star_allsame'] = all_same
    if all_same:
        if not np.allclose(df['mu_star_flowline_avg'],
                           df['mu_star_glacierwide'],
                           atol=1e-3):
            raise MassBalanceCalibrationError('Unexpected difference between '
                                              'glacier wide mu* and the '
                                              'flowlines mu*.')
    # Write
    gdir.write_json(df, 'local_mustar')
```
As a final step, the following diagnostic parameters are written to the `local_mustar.json`:

- `mu_star_per_flowline` is an array containing the temperature sensitivity $\mu^*$ for each flow line.
- `mu_star_flowline_avg` is the average $\overline{\mu^*}$ over all flow lines, weighted per width.
- `mu_star_allsame` is `True` if $\mu^*$ of all flow lines lie within a range of $0.001$ and `False` otherwise.

If $\mu^*$ of all flow lines are equal but differ significantly from the glacier wide $\mu^*$ an error is raised.



### `_recursive_mu_star_calibration(gdir, ...)`

This function is called by the `mu_star_calibration(gdir)` task and (as the name may suggest) recursively calibrates the $\mu^*$ for each flow line, accounting for possibly arising negative fluxes.

```python
def _recursive_mu_star_calibration(gdir, fls, t_star, first_call=True):
```

The function works on the GlacierDirectory `gdir`. In addition the glacier flow lines `fls` and $t^*$ are needed. The flag `first_call` is set to `True` by default, disabling certain computation on the following recursive calls.

```python
    # Do we have a calving glacier? This is only for the first call!
    # The calving mass-balance is distributed over the valid tributaries of the
    # main line, i.e. bad tributaries are not considered for calving
    cmb = calving_mb(gdir) if first_call else 0.
```

On the first call the function checks for a possible calving mass balance `cmb`.

```python
    # Climate period
    mu_hp = int(cfg.PARAMS['mu_star_halfperiod'])
    yr_range = [t_star - mu_hp, t_star + mu_hp]

    # Get the corresponding mu
    heights = np.array([])
    widths = np.array([])
    for fl in fls:
        heights = np.append(heights, fl.surface_h)
        widths = np.append(widths, fl.widths)
        
    _, temp, prcp = mb_yearly_climate_on_height(gdir, heights,
                                                year_range=yr_range,
                                                flatten=False)
```

As before, the climate period is defined as 31 years centred around $t^*$. The elevation levels and corresponding glacier widths of each flow line are stored in the `heights` and `widths` container, respectively. The `mb_yearly_climate_on_height(gdir)` returns the sum of the driving meteorological parameters over a year on the given heights, as discussed above.

```python
    try:
        mu_star = optimization.brentq(_mu_star_per_minimization,
                                      cfg.PARAMS['min_mu_star'],
                                      cfg.PARAMS['max_mu_star'],
                                      args=(fls, cmb, temp, prcp, widths),
                                      xtol=1e-5)
    except ValueError:
        # TODO: add "f(a) and f(b) must have different signs" check
        raise MassBalanceCalibrationError('{} has mu which exceeds the '
                                          'specified min and max '
                                          'boundaries.'.format(gdir.rgi_id))

    if not np.isfinite(mu_star):
        raise MassBalanceCalibrationError('{} '.format(gdir.rgi_id) +
                                          'has a non finite mu.')


```

The `optimization.brentq` function tries to find a root of the given function using Brent’s method (whatever that is, [see docs for detail](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.brentq.html)). Here it is used to optimise the $\mu^*$ parameter so that the mass balance over the climate period centred around $t^*$ is zero, with a tolerance of $10^{-5}$. Thereby the `_mu_star_per_minimizatuion(...)` returns a value for said mass balance over the climate period, see below for details. The surrounding `try catch` block raises an error if the optimised $\mu^*$ lies outside the bounds specified in the `PARAMS.cfg` file. In case of a non finite result a fitting error is raised.

```python
    # Reset flux
    for fl in fls:
        fl.flux = np.zeros(len(fl.surface_h))

    # Flowlines in order to be sure - start with first guess mu*
    for fl in fls:
        y, t, p = mb_yearly_climate_on_height(gdir, fl.surface_h,
                                              year_range=yr_range,
                                              flatten=False)
        mu = fl.mu_star if fl.mu_star_is_valid else mu_star
        fl.set_apparent_mb(np.mean(p, axis=1) - mu*np.mean(t, axis=1),
                           mu_star=mu)
```

Next step is to calculate the apparent mass balance and the corresponding flux for each flow line. Before doing so the flux of each flow line must be reset to zero. The mass balance is given as by equation (see theoretical background) as the difference between average solid precipitation (mass gain) and the average positive temperature times the temperature sensitivity $\mu^*$ (mass loss).

```python
    # Sometimes, low lying tributaries have a non-physically consistent
    # Mass-balance. These tributaries wouldn't exist with a single
    # glacier-wide mu*, and therefore need a specific calibration.
    # All other mus may be affected
    if cfg.PARAMS['correct_for_neg_flux']:
        if np.any([fl.flux_needs_correction for fl in fls]):

            # We start with the highest Strahler number that needs correction
            not_ok = np.array([fl.flux_needs_correction for fl in fls])
            fl = np.array(fls)[not_ok][-1]

            # And we take all its tributaries
            inflows = centerlines.line_inflows(fl)

            # We find a new mu for these in a recursive call
            # TODO: this is where a flux kwarg can passed to tributaries
            _recursive_mu_star_calibration(gdir, inflows, t_star,
                                           first_call=False)

            # At this stage we should be ok
            assert np.all([~ fl.flux_needs_correction for fl in inflows])
            for fl in inflows:
                fl.mu_star_is_valid = True

            # After the above are OK we have to recalibrate all below
            _recursive_mu_star_calibration(gdir, fls, t_star,
                                           first_call=first_call)
            
    # At this stage we are good
    for fl in fls:
        fl.mu_star_is_valid = True
```

Here the recursion kicks in, if the `correct_for_neg_flux` in the `PARAMS.cfg` is set. This correction is necessary since, sometimes low lying (low in altitude) tributaries have a non-physically consistent mass balance. These tributaries wouldn't exist with a single glacier wide $\mu^*$, and therefore need a specific calibration. The $\mu^*$ of all other flow lines may be affected.

Flow lines with negative fluxes in the last two grid points need correction, starting with the highest Strahler number (i.e. the 'most parent' flow line). All its tributaries are passed to the `_recursive_mu_star_calibration(gdir, ...)`, which calibrates a new $\mu^*$ for the inflows. Arising negative fluxes in tributaries during this secondary (tertiary, ...) calibration are caught by the same mechanism. When all current (i.e. erroneous) flow lines are calibrated correctly, the calibration process is repeated with all underlying flow lines (i.e. more 'parent' flow lines).

If the recursive process comes to an end, the calibration process is done.



### `_mu_star_per_minimization(...)`

This function is used by the `_recursive_mu_star_calibration(...)` task via a root finding method.

```python
def _mu_star_per_minimization(x, fls, cmb, temp, prcp, widths):
```

The parameters of this function are defined to be used with a root finding method, whereby `x` is the guessed value for the temperature sensitivity $\mu^*$. In addition the glacier flow lines `fls`, calving mass balance `cmb`, yearly positive temperature and solid precipitation on the flow line elevation levels `temp` and `prcp`, as well as the corresponding glacier widths `widths` are needed.

```python
    # Get the corresponding mu
    mus = np.array([])
    for fl in fls:
        mu = fl.mu_star if fl.mu_star_is_valid else x
        mus = np.append(mus, np.ones(fl.nx) * mu)
```

For each glacier flow line the $\mu^*$ is set to the optimisation parameter `x` unless the `mu_star_is_valid` flag is set. In this case the value stored in the flow line object `fl.mu_star` is used. The value is stored as an array with the same number of elements as the elevation levels of the corresponding flow line.

```python
    # TODO: possible optimisation here
    out = np.average(prcp - mus[:, np.newaxis] * temp, axis=0, weights=widths)
```

The mass balance is computed for every elevation level over the climate period for each flow line. Afterwards it is averaged weighted on glacier width and stored in the `out` parameter. Given the definition of $t^*$, the sum over the climate period should be zero (therefore the root finding function).

```python
    return np.mean(out - cmb)
```

Before returning the average over all flow lines is computed and it is accounted for a possible calving mass balance.



### `Centerline.set_apparent_mb(mb)`

This method is used during the $\mu^*$ calibration by the `_recursive_mu_star_calibration(...)` task.

```python
def set_apparent_mb(self, mb, mu_star=None):
    """Set the apparent mb and flux for the flowline.

    MB is expected in kg m-2 yr-1 (= mm w.e. yr-1)

    This should happen in line order, otherwise it will be wrong.

    Parameters
    ----------
    mu_star : float
        if appropriate, the mu_star associated with this apparent mb
    """
```

This is a method of the `Centerline` class, setting the apparent mass balance `mb` (in $\mathrm{kg\ m^{-2}\ yr^{-1}}$ ) and ice flux for the calling flow line. It is possible to supply the corresponding $\mu^*$.

```python
	self.apparent_mb = mb
    self.mu_star = mu_star
```

The given mass balance and $\mu^*$ are just stored as flow line parameters.

```python
    # Add MB to current flux and sum
    # no more changes should happen after that
    flux_needs_correction = False
    flux = np.cumsum(self.flux + mb * self.widths * self.dx)
    
	# We filter only lines with two negative grid points, the
    # rest we can cope with
    if flux[-2] < 0:
        flux_needs_correction = True
        
	self.flux = flux
    self.flux_needs_correction = flux_needs_correction
```

The flux for every flow line element is calculated from the mass balance and summed up. If the last two grid points show a negative flux some correction is needed. Values are stored as flow line attributes.

```python
    # Add to outflow. That's why it should happen in order
    if self.flows_to is not None:
        n = len(self.flows_to.line.coords)
        ide = self.flows_to_indice
        if n >= 9:
            gk = GAUSSIAN_KERNEL[9]
            self.flows_to.flux[ide-4:ide+5] += gk * flux[-1]
        elif n >= 7:
            gk = GAUSSIAN_KERNEL[7]
            self.flows_to.flux[ide-3:ide+4] += gk * flux[-1]
        elif n >= 5:
            gk = GAUSSIAN_KERNEL[5]
            self.flows_to.flux[ide-2:ide+3] += gk * flux[-1]
```

If the given `Centerline` flows into another one, its outflowing mass must be *transported* to the parent flow line. This is done by adding the outflow, i.e. the flux of the last grid points `flux[-1]`, to the flux of the parent flow line. The flux is distributed in over several grid points centred at the intersection point using a Gaussian distribution. The width of the Gaussian distribution depends on the lengths of the parent flow line. This is why the apparent mass balance computation must be done in order, starting with tributaries and ending with the main flow line.



## The `PARAMS.cfg` file

Hereafter copied are all above used parameters from the default `PARAMS.cfg` file with the explaining comments.

```python
# specify here the start and end year where oggm will searh for tstar
# candidates (note that the window will be reduced by mu_star_halfperiod on
# each side of the window). Set to 0, 0 for the default (the entire available
# data space)
tstar_search_window = 0, 0
mu_star_halfperiod = 15

# which temperature gradient? if false, use temp_default_gradient. If true,
# compute by regression of the 9 surrounding grid points (not recommended)
temp_use_local_gradient = False
temp_default_gradient = -0.0065
# the linear regression can lead to quite strange results... this helps
# you to clip them to more realistic values:
temp_local_gradient_bounds = -0.009, -0.003
# other parameters
temp_all_solid = 0.
temp_all_liq = 2.
temp_melt = -1.
# precipitation correction: set to a float for a constant scaling factor
prcp_scaling_factor = 2.5

# Should we use the default, pre-calibrated reference tstars or are we
# running the calibration ourselves? The default should be False, which
# raises a warning when trying to calibrate.
run_mb_calibration = False
# Bounds on mu*
# Values out of these limits are considered bad and will lead to an error
min_mu_star = 0.
max_mu_star = 10000.
```

