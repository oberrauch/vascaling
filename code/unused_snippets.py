def get_climatological_temp_prcp(gdir, time_range=None, year_range=None):
    """ Get the monthly climatological values around t* of solid precipitation and temperature.
    ATTENTION: Works, but should not be used for mu* calibration. TODO: Will/Should be deleted.

    :param gdir:
    :param min_hgt:
    :param max_hgt:
    :param time_range:
    :param year_range:
    :return:
    """
    if year_range is not None:
        sm = cfg.PARAMS['hydro_month_' + gdir.hemisphere]
        em = sm - 1 if (sm > 1) else 12
        t0 = datetime.datetime(year_range[0]-1, sm, 1)
        t1 = datetime.datetime(year_range[1], em, 1)
        return get_climatological_temp_prcp(gdir, time_range=[t0, t1])

    # Parameters
    temp_all_solid = cfg.PARAMS['temp_all_solid']
    prcp_fac = cfg.PARAMS['prcp_scaling_factor']
    default_grad = cfg.PARAMS['temp_default_gradient']
    prcp_grad = 3e-4
    g_minmax = cfg.PARAMS['temp_local_gradient_bounds']

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
                raise climate.MassBalanceCalibrationError('time_range[0] not found in '
                                                          'file')
            p1 = np.where(time == time_range[1])[0]
            try:
                p1 = p1[0]
            except IndexError:
                raise climate.MassBalanceCalibrationError('time_range[1] not found in '
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

    # The following is my code. So abandon all hope, you who enter here.

    # get relevant elevation information
    fpath = gdir.get_filepath('gridded_data')
    with ncDataset(fpath) as nc:
        mask = nc.variables['glacier_mask'][:]
        topo = nc.variables['topo'][:]
    min_elev = np.min(topo[np.where(mask == 1)])
    max_elev = np.max(topo[np.where(mask == 1)])

    # get temperature at glacier terminus
    temp_terminus = _compute_temp_terminus(itemp, igrad, ref_hgt, min_elev)
    # get solid precipitation
    prcp_solid = _compute_solid_prcp(iprcp, prcp_fac, ref_hgt, min_elev, max_elev,
                                     temp_terminus, temp_all_solid, igrad, prcp_grad)

    # compute climatological values
    months = np.arange(1, 13)
    temp_terminus_clim = months * 0
    prcp_solid_clim = months * 0
    for month in months:
        # get location
        pos = [t.month == month for t in time]
        # compute average over all months
        prcp_solid_clim[month - 1] = np.mean(prcp_solid[pos])
        temp_terminus_clim[month - 1] = np.mean(temp_terminus[pos])

    return months, temp_terminus_clim, prcp_solid_clim

class OriginalMassBalance(MassBalanceModel):
    """Original mass balance model from Ben's paper.
    As as starting point, I just copied the PastMassBalance model.
    """

    def __init__(self, gdir, mu_star=None, bias=None,
                 filename='climate_monthly', input_filesuffix='',
                 repeat=False, ys=None, ye=None, check_calib_params=True):
        """Initialize. Exactly the same as in the PastMassBalance model (for now).

        Parameters
        ----------
        gdir : GlacierDirectory
            the glacier directory
        mu_star : float, optional
            set to the alternative value of mu* you want to use
            (the default is to use the calibrated value).
        bias : float, optional
            set to the alternative value of the calibration bias [mm we yr-1]
            you want to use (the default is to use the calibrated value)
            Note that this bias is *substracted* from the computed MB. Indeed:
            BIAS = MODEL_MB - REFERENCE_MB.
        filename : str, optional
            set to a different BASENAME if you want to use alternative climate
            data.
        input_filesuffix : str
            the file suffix of the input climate file
        repeat : bool
            Whether the climate period given by [ys, ye] should be repeated
            indefinitely in a circular way
        ys : int
            The start of the climate period where the MB model is valid
            (default: the period with available data)
        ye : int
            The end of the climate period where the MB model is valid
            (default: the period with available data)
        check_calib_params : bool
            OGGM will try hard not to use wrongly calibrated mu* by checking
            the parameters used during calibration and the ones you are
            using at run time. If they don't match, it will raise an error.
            Set to False to suppress this check.

        Attributes
        ----------
        temp_bias : float, default 0
            Add a temperature bias to the time series
        prcp_bias : float, default 1
            Precipitation factor to the time series (called bias for
            consistency with `temp_bias`)
        """

        # call the initialisator from the parent class
        super(OriginalMassBalance, self).__init__()
        # specify valid elevation bounds in meters
        self.valid_bounds = [-1e4, 2e4]

        # read mu* from parameter file if not given
        if mu_star is None:
            df = gdir.read_json('ben_params')
            mu_star = df['mu_star']

        # read bias from parameter file if not given
        if bias is None:
            if cfg.PARAMS['use_bias_for_run']:
                df = gdir.read_json('ben_params')
                bias = df['bias']
            else:
                bias = 0.

        # set mass balance attributes
        self.mu_star = mu_star
        self.bias = bias

        # Parameters
        self.t_solid = cfg.PARAMS['temp_all_solid']
        self.t_liq = cfg.PARAMS['temp_all_liq']
        self.t_melt = cfg.PARAMS['temp_melt']
        prcp_fac = cfg.PARAMS['prcp_scaling_factor']
        default_grad = cfg.PARAMS['temp_default_gradient']

        # Check the climate related params to the GlacierDir to make sure
        if check_calib_params:
            mb_calib = gdir.read_pickle('climate_info')['mb_calib_params']
            for k, v in mb_calib.items():
                if v != cfg.PARAMS[k]:
                    raise RuntimeError('You seem to use different mass-'
                                       'balance parameters than used for the '
                                       'calibration. '
                                       'Set `check_calib_params=False` '
                                       'to ignore this warning.')

        # Public attributes
        self.temp_bias = 0.
        self.prcp_bias = 1.
        self.repeat = repeat

        # Read file
        fpath = gdir.get_filepath(filename, filesuffix=input_filesuffix)
        with ncDataset(fpath, mode='r') as nc:
            # time
            time = nc.variables['time']
            time = netCDF4.num2date(time[:], time.units)
            # check if climate data includes all 12 month for each year
            ny, r = divmod(len(time), 12)
            if r != 0:
                raise ValueError('Climate data should be N full years')
            # This is where we switch to hydro float year format
            # Last year gives the tone of the hydro year
            self.years = np.repeat(np.arange(time[-1].year-ny+1,
                                             time[-1].year+1), 12)
            self.months = np.tile(np.arange(1, 13), ny)
            # Read time series
            self.temp = nc.variables['temp'][:]
            self.prcp = nc.variables['prcp'][:] * prcp_fac
            if 'gradient' in nc.variables:
                grad = nc.variables['gradient'][:]
                # Security for stuff that can happen with local gradients
                g_minmax = cfg.PARAMS['temp_local_gradient_bounds']
                grad = np.where(~np.isfinite(grad), default_grad, grad)
                grad = np.clip(grad, g_minmax[0], g_minmax[1])
            else:
                grad = self.prcp * 0 + default_grad
            self.grad = grad
            self.grad_prcp = 3e-4
            self.ref_hgt = nc.ref_hgt
            self.ys = self.years[0] if ys is None else ys
            self.ye = self.years[-1] if ye is None else ye

    def get_monthly_climate(self, heights, year=None):
        """ Gather monthly climate information,
        i.e. compute terminus temperature and solid precipitation.

        Note that prcp is corrected with the precipitation factor (already)
        and that all other model biases (temp and prcp) are applied.

        Returns
        -------
        (temp_terminus, prcpsol)
        """

        # get year and month from decimal year (fraction)
        y, m = floatyear_to_date(year)
        # deal with years outside the period with climate data
        if self.repeat:
            y = self.ys + (y - self.ys) % (self.ye - self.ys + 1)
        if y < self.ys or y > self.ye:
            raise ValueError('year {} out of the valid time bounds: '
                             '[{}, {}]'.format(y, self.ys, self.ye))

        # get position in time series corresponding
        # to given year and month combination
        pok = np.where((self.years == y) & (self.months == m))[0][0]

        # read temperature and precipitation time series
        # and apply possible bias
        itemp = self.temp[pok] + self.temp_bias
        itemp_anom = 0
        iprcp = self.prcp[pok] * self.prcp_bias
        iprcp_anom = 0
        # read temperature gradient
        igrad = self.grad[pok]

        # @WHERE_TO_START

        # compute the temperature at the terminus and the highest point
        temp_terminus = itemp + igrad * (min(heights) - self.ref_hgt) + itemp_anom
        temp_zmax = temp_terminus + igrad * (max(heights) - min(heights))

        # compute solid precipitation fraction
        f_solid = 1 + ((temp_terminus - self.t_solid) / (igrad * (max(heights) - min(heights))))
        if temp_terminus <= self.t_solid:
            f_solid = 1
        elif temp_zmax >= self.t_solid:
            f_solid = 0

        # compute solid precipitation
        prcpsol = iprcp + iprcp_anom
        prcpsol *= (1 + self.grad_prcp * (np.mean([heights[0], heights[-1]]) - self.ref_hgt))
        prcpsol *= f_solid

        return temp_terminus, prcpsol

    def get_monthly_mb(self, heights, year=None, fl_id=None):
        # get terminus temperature and solid precipitation for given year/month
        temp_terminus, prcpsol = self.get_monthly_climate(heights, year=year)
        # compute monthly mass balance
        mb_month = prcpsol - self.mu_star * max((temp_terminus - self.t_melt), 0)
        # apply bias (scaled per month)
        mb_month -= self.bias * SEC_IN_MONTH / SEC_IN_YEAR
        # convert into SI units
        return mb_month / SEC_IN_MONTH / self.rho

    def get_annual_mb(self, heights, year=None, fl_id=None):
        # create float dates for all month in given year
        months = np.arange(year, year+1, 1/12)
        # sum monthly mb over the given year
        mb_year = 0
        for month in months:
            mb_year += self.get_monthly_mb(heights, month, fl_id)
        return mb_year * SEC_IN_MONTH / SEC_IN_YEAR

    def run_and_store(self, year_end, reset=False):
        """ Runs the model till the specified year.
        Returns all geometric parameters (i.e. lenght, area, volume and
        terminus elevation) for each modelled year, as well as the years as
        index. TODO: run_and_store() or similar, but for now ok?!

        :param year_end: (float) end of modeling period
        :param reset: (bool, optional) If `True`, the model will start from
            `year_0`, otherwise from its current position in time (default).
        :return:
        """

        # reset parameters to starting values
        if reset:
            self.reset()

        # check validity of end year
        if year_end < self.year:
            raise ValueError('Cannot run until {}, already at year {}'.format(
                year_end, self.year))

        # create empty containers
        years = np.arange(self.year, year_end + 1)
        year0 = np.empty(years.size)
        length_m = np.empty(years.size)
        length_m_0 = np.empty(years.size)
        dL = np.empty(years.size)

        area_m2 = np.empty(years.size)
        area_m2_0 = np.empty(years.size)
        dA = np.empty(years.size)

        volume_m3 = np.empty(years.size)
        volume_m3_0 = np.empty(years.size)
        dV = np.empty(years.size)

        tau_l = np.empty(years.size)
        tau_a = np.empty(years.size)

        ca = np.empty(years.size)
        gamma = np.empty(years.size)
        cl = np.empty(years.size)
        ql = np.empty(years.size)

        min_hgt = np.empty(years.size)
        min_hgt_0 = np.empty(years.size)
        max_hgt = np.empty(years.size)
        rho = np.empty(years.size)
        spec_mb = np.empty(years.size)

        # iterate over all years
        for i, year in enumerate(years):
            if i != 0:
                # run model for one year
                self.step()

            # store metrics
            year0[i] = self.year_0
            
            length_m[i] = self.length_m
            length_m_0[i] = self.length_m_0
            dL[i] = self.dL

            area_m2[i] = self.area_m2
            area_m2_0[i] = self.area_m2_0
            dA[i] = self.dA

            volume_m3[i] = self.volume_m3
            volume_m3_0[i] = self.volume_m3_0
            dV[i] = self.dV

            tau_l[i] = self.tau_l
            tau_a[i] = self.tau_a

            ca[i] = self.ca
            gamma[i] = self.gamma
            cl[i] = self.cl
            ql[i] = self.ql

            min_hgt[i] = self.min_hgt
            min_hgt_0[i] = self.min_hgt_0
            max_hgt[i] = self.max_hgt
            rho[i] = self.rho
            spec_mb[i] = self.spec_mb

        # create DataFrame
        columns = ['year0', 'length_m', 'length_m_0', 'dL', 'area_m2',
                   'area_m2_0', 'dA', 'volume_m3', 'volume_m3_0', 'dV',
                   'tau_l', 'tau_a', 'ca', 'gamma', 'cl', 'ql', 'min_hgt',
                   'min_hgt_0', 'max_hgt', 'rho', 'spec_mb']
        df = pd.DataFrame(np.array([year0, length_m, length_m_0, dL, area_m2,
                                    area_m2_0, dA, volume_m3, volume_m3_0, dV,
                                    tau_l, tau_a, ca, gamma, cl, ql, min_hgt,
                                    min_hgt_0, max_hgt, rho, spec_mb]).T,
                          index=years, columns=columns)
        # return metrics
        return df
