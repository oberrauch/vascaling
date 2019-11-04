"""
    Hereafter I'll try to set up one run with my new model.
    This includes:
    a) initialisation and calibration
    b) the mass balance model
    c) the 'dynamic' model

    Date: 18.01.2019

    Update on January, 16: Run automated comparison between both models,
    including initialization, computation, plots, ...
    -> New features to implement thereby:
        [ ] only RGI ID as input
        [x] use 'vas_ref_tstars.csv'
        [ ] compute starting area

"""
# import externals libs
import os
import shutil
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt


# import the needed OGGM modules
import oggm
from oggm import cfg, utils
from oggm.core import (gis, climate, centerlines, flowline,
                       inversion, massbalance, vascaling)


def compare(rgi_id, glacier_name):
    """

    :param rgi_id:
    :param glacier_name:
    :return:
    """

    # ---------------------
    #  PREPROCESSING TASKS
    # ---------------------
    # create test directory
    wdir = os.path.join(os.path.abspath('.'), 'comparison_wdir')
    if not os.path.exists(wdir):
        os.makedirs(wdir)
    shutil.rmtree(wdir)
    os.makedirs(wdir)

    # load default parameter file
    cfg.initialize()

    # RGI entity
    # get/downlaod the rgi entity including the outline shapefile
    rgi_df = utils.get_rgi_glacier_entities([rgi_id])
    # set name, since not delivered with RGI
    if rgi_df.loc[int(rgi_id[-5:])-1, 'Name'] is None:
        rgi_df.loc[int(rgi_id[-5:])-1, 'Name'] = glacier_name

    # select single entry
    rgi_entity = rgi_df.iloc[0]

    # GlacierDirectory
    # specify the working directory and define the glacier directory
    cfg.PATHS['working_dir'] = wdir
    gdir = oggm.GlacierDirectory(rgi_entity)

    # DEM and GIS tasks
    # get the path to the DEM file (will download if necessary)
    dem = utils.get_topo_file(gdir.cenlon, gdir.cenlat)
    # set path in config file
    cfg.PATHS['dem_file'] = dem[0][0]
    cfg.PARAMS['border'] = 10
    cfg.PARAMS['use_intersects'] = False
    # run GIS tasks
    gis.define_glacier_region(gdir, entity=rgi_entity)
    gis.glacier_masks(gdir)

    # Climate data
    # using HistAlp
    cfg.PARAMS['baseline_climate'] = 'HISTALP'
    # climate records before 1850 are hardly reliable, which is not so drastic for
    # qualitative experiments (could be driven with random climate anyway)
    # cfg.PARAMS['baseline_y0'] = 1850
    # change hyper parameters for HistAlp
    cfg.PARAMS['prcp_scaling_factor'] = 1.75
    cfg.PARAMS['temp_melt'] = -1.75
    # run climate task
    climate.process_histalp_data(gdir)

    # run center line preprocessing tasks
    centerlines.compute_centerlines(gdir)
    centerlines.initialize_flowlines(gdir)
    centerlines.compute_downstream_line(gdir)
    centerlines.compute_downstream_bedshape(gdir)
    centerlines.catchment_area(gdir)
    centerlines.catchment_intersections(gdir)
    centerlines.catchment_width_geom(gdir)
    centerlines.catchment_width_correction(gdir)

    # --------------------
    #  SCALING MODEL
    # --------------------

    # compute local t* and the corresponding mu*
    vascaling.local_t_star(gdir)

    # instance the mass balance models
    vas_mb_mod = vascaling.VAScalingMassBalance(gdir)

    # get reference area
    a0 = gdir.rgi_area_m2
    # get reference year
    y0 = gdir.read_pickle('climate_info')['baseline_hydro_yr_0']
    y1 = gdir.read_pickle('climate_info')['baseline_hydro_yr_1']
    # get min and max glacier surface elevation
    h0, h1 = vascaling.get_min_max_elevation(gdir)

    # instance VAS model
    vas_model = vascaling.VAScalingModel(year_0=y0, area_m2_0=a0,
                                         min_hgt=h0, max_hgt=h1,
                                         mb_model=vas_mb_mod)
    # run model over all HistAlp climate period
    vas_df = vas_model.run_and_store(y1, reset=True)
    # get relevant parameters
    years_vas = vas_df.index.values
    length_m_vas = vas_df.length_m.values
    area_m2_vas = vas_df.area_m2.values
    volume_m3_vas = vas_df.volume_m3.values

    # ------
    #  OGGM
    # ------

    # compute local t* and the corresponding mu*
    climate.local_t_star(gdir)
    climate.mu_star_calibration(gdir)

    # instance the mass balance models
    mb_mod = massbalance.PastMassBalance(gdir)

    # run inversion tasks
    inversion.prepare_for_inversion(gdir)
    inversion.mass_conservation_inversion(gdir)
    inversion.filter_inversion_output(gdir)

    # initialize present time glacier
    flowline.init_present_time_glacier(gdir)

    # instance flowline model
    fls = gdir.read_pickle('model_flowlines')
    y0 = gdir.read_pickle('climate_info')['baseline_hydro_yr_0']
    y1 = gdir.read_pickle('climate_info')['baseline_hydro_yr_1']
    fl_mod = flowline.FluxBasedModel(flowlines=fls, mb_model=mb_mod, y0=y0)

    # run model and store output as xarray data set
    _, oggm_ds = fl_mod.run_until_and_store(y1)

    years_oggm = oggm_ds.hydro_year.values
    # annual index must be equal
    np.testing.assert_array_equal(years_oggm, years_vas)
    length_m_oggm = oggm_ds.length_m.values
    area_m2_oggm = oggm_ds.area_m2.values
    volume_m3_oggm = oggm_ds.volume_m3.values

    # define column names for DataFrame
    names = ['length_vas', 'length_oggm',
             'area_vas', 'area_oggm',
             'volume_vas', 'volume_oggm']
    # combine glacier geometries into DataFrame
    df = pd.DataFrame(np.array([length_m_vas, length_m_oggm,
                                area_m2_vas, area_m2_oggm,
                                volume_m3_vas, volume_m3_oggm]).T,
                      index=years_vas, columns=names)
    # save to file
    store = True
    if store:
        # define path and file names
        folder = '/Users/oberrauch/work/master/data/'
        df.to_csv(folder+'run_comparison.csv')

    def plot_both(vas_df, oggm_df, ref=None, correct_bias=False,
                  title='', ylabel='', file_path='', exp=0):
        """ Plot geometric parameters of both models.
        If a `file_path` is given, the figure will be saved.

        :param vas_df: (pandas.Series) geometric glacier parameter of the VAS model
        :param oggm_df: (pandas.Series) geometric glacier parameter of the OGGM
        :param ref: (pandas.Series) measured glacier parameter, optional
        :param title: (string) figure title, optional
        :param ylabel: (string) label for y-axis, optional
        :param file_path: (string) where to store the figure, optional
        :param exp: (int) exponent for labels in scientific notation, optional
        """
        beamer = True
        if beamer:
            mpl.rc('axes', titlesize=18)
            mpl.rc('axes', labelsize=14)
            mpl.rc('xtick', labelsize=14)
            mpl.rc('ytick', labelsize=14)
            mpl.rc('legend', fontsize=10)
        # create figure and first axes
        fig = plt.figure(figsize=[6, 4])
        ax = fig.add_axes([0.15, 0.1, 0.8, 0.8])

        # define colors
        c1 = 'C0'
        c2 = 'C1'
        c3 = 'C3'
        # plot vas and OGGM parameters
        ax.plot(oggm_df.index, oggm_df.values, c=c2, label='OGGM')
        ax.plot(vas_df.index, vas_df.values, c=c1, label='VAS')
        if ref:
            # plot reference parameter if given
            ax.plot(ref.index, ref.values, c=c3, label='measurements')
        if correct_bias:
            # plot bias corrected vas
            df_ = pd.DataFrame([oggm_df, vas_df]).T
            bias = vas_df.values - df_.mean().diff().iloc[1]
            ax.plot(vas_df.index, bias, c=c1, ls='--',
                    label='VAS, bias corrected')
            # add RMSD as text
            ax.text(0.05, 0.05,
                    'RMSD: {:.1e}'.format(utils.rmsd(oggm_df, bias)),
                    transform=plt.gca().transAxes)

        # add correlation coefficient as text
        ax.text(0.05, 0.11, 'Corr. Coef.: {:.2f}'.format(
            utils.corrcoef(oggm_df, vas_df)),
                    transform=plt.gca().transAxes)

        # add title, labels, legend
        ax.set_title(title)
        ax.set_ylabel(ylabel)
        ax.legend()

        import matplotlib.ticker

        class OOMFormatter(matplotlib.ticker.ScalarFormatter):
            def __init__(self, order=0, fformat="%1.1f", offset=False, mathText=False):
                self.oom = order
                self.fformat = fformat
                matplotlib.ticker.ScalarFormatter.__init__(self, useOffset=offset, useMathText=mathText)

            def _set_orderOfMagnitude(self, nothing):
                self.orderOfMagnitude = self.oom

            def _set_format(self, vmin, vmax):
                self.format = self.fformat
                if self._useMathText:
                    self.format = '$%s$' % matplotlib.ticker._mathdefault(self.format)

        # use scientific notation with fixed exponent according
        ax.yaxis.set_major_formatter(OOMFormatter(exp, "%1.2f"))

        # store to file
        if file_path:
            plt.savefig(file_path, bbox_inches='tight',
                        format=file_path.split('.')[-1])

    # specify plot directory
    folder = '/Users/oberrauch/work/master/plots/'

    # plot length
    plot_both(df.length_vas, df.length_oggm, correct_bias=True,
              title='Glacier length - {}'.format(glacier_name),
              ylabel=r'Length [m]',
              file_path=os.path.join(folder, '{}_length.pdf'.format(rgi_id)),
              exp=3)
    # plot area
    plot_both(df.area_vas, df.area_oggm, correct_bias=True,
              title='Surface area - {}'.format(glacier_name),
              ylabel=r'Area [m$^2$]',
              file_path=os.path.join(folder, '{}_area.pdf'.format(rgi_id)),
              exp=6)
    # plot volume
    plot_both(df.volume_vas, df.volume_oggm, correct_bias=True,
              title='Glacier volume - {}'.format(glacier_name),
              ylabel=r'Volume [m$^3$]',
              file_path=os.path.join(folder, '{}_volume.pdf'.format(rgi_id)),
              exp=9)


if __name__ == '__main__':
    # choose glacier
    rgi_id = 'RGI60-11.00897'
    glacier_name = 'Hintereisferner'
    compare(rgi_id, glacier_name)
