# import section
import os
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

from plots.master_colors import vas_cycle, fl_cycle
# from master_colors import vas_cycle, fl_cycle

import logging

logging.basicConfig(format='%(asctime)s: %(name)s: %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    level=logging.INFO)
log = logging.getLogger('plot timeseries')


def read_dataset(path):
    # read dataset
    ds = xr.load_dataset(path)
    # sort by temperature bias
    ds = ds.sortby('temp_bias')
    # cast normalized dimension from int to bool
    ds['normalized'] = [bool(norm) for norm in ds.normalized]

    return ds


# ---------------------------------------------
# Defining plotting functions for universal use
# ---------------------------------------------


def plot_time_series(ds, var, x='time',
                     title='', xlabel='Years of model evolution', ylabel='',
                     labels=None, legend_title='', legend_loc=0,
                     color_cycle=None, xlim=None, path=None):
    # plot relative volume change
    fig, ax = plt.subplots(1, 1, figsize=[5, 4])

    # set color cycle
    if color_cycle is not None:
        ax.set_prop_cycle('color', color_cycle)

    # plot
    handles = ds[var].plot.line(x=x, add_legend=False)

    # add legend
    if labels:
        legend = ax.legend(handles, labels, title=legend_title, loc=legend_loc)

    # title, axis labels, ...
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    if not ylabel:
        if ds.normalized:
            # add ylabel
            ax.set_ylabel('Relative {}'.format(var))
            # aux line
            ax.axhline(1, lw=0.8, ls=':', c='k')
        else:
            # add ylabel
            unit = 'm' if var == 'length' else 'm$^2$' if var == 'area' else 'm$^3$'
            ax.set_ylabel('Glacier {} [{}]'.format(var, unit))
            # aux line
            ax.axhline(ds[var].isel(time=0).mean(), lw=0.8, ls=':', c='k')
    else:
        ax.set_ylabel(ylabel)
    ax.grid()
    # set axes limits
    if xlim:
        ax.set_xlim(xlim)
    else:
        ax.set_xlim([ds[x].min(), ds[x].max()])

    if path:
        fig.savefig(path, bbox_inches='tight')


def plot_eq_time_series(ds, var, title='', suptitle='', xlim=None):
    # plot relative volume change
    fig, [ax0, ax1] = plt.subplots(1, 2, figsize=[10, 5])

    # flowline model
    ds.sel(model='fl')[var].plot(hue='temp_bias', ax=ax0, add_legend=False,
                                 color='lightgray', lw=0.5)
    # vas model
    ax0.set_prop_cycle('color', vas_cycle)
    handles_vas = ds.sel(model='vas')[var].plot(hue='temp_bias', ax=ax0,
                                                add_legend=False, lw=2)
    labels_vas = ['{:+.1f} °C'.format(bias) for bias in ds.temp_bias.values]

    # set axes limits
    if xlim:
        ax0.set_xlim(xlim)
    else:
        ax0.set_xlim([ds.time.min(), ds.time.max()])
    ylim = ax0.get_ylim()
    ax0.set_ylim([ylim[0] * 0.5, ylim[1]])

    # title, labels, legend
    ax0.set_title('Volume/area scaling model')
    ax0.set_xlabel('Years of model evolution')
    ax0.legend(handles_vas, labels_vas, title='Temperature bias',
               bbox_to_anchor=(0.5, 0), loc=8, ncol=3)
    if ds.normalized:
        # add ylabel
        ax0.set_ylabel('Relative {}'.format(var))
        # aux line
        ax0.axhline(1, lw=0.8, ls=':', c='k')
    else:
        # add ylabel
        unit = 'm' if var == 'length' else 'm$^2$' if var == 'area' else 'm$^3$'
        ax0.set_ylabel('Glacier {} [{}]'.format(var, unit))
        # aux line
        ax0.axhline(ds.sel(model='vas')[var].isel(time=0).mean(), lw=0.8,
                    ls=':', c='k')

    ax0.grid()

    # vas model
    ds.sel(model='vas')[var].plot(hue='temp_bias', ax=ax1, add_legend=False,
                                  color='lightgray', lw=0.5)
    # flowline model
    ax1.set_prop_cycle('color', fl_cycle)
    handles_fl = ds.sel(model='fl')[var].plot(hue='temp_bias', ax=ax1,
                                              add_legend=False)
    labels_fl = ['{:+.1f} °C'.format(bias) for bias in ds.temp_bias.values]

    # set axes limits
    if xlim:
        ax1.set_xlim(xlim)
    else:
        ax1.set_xlim([ds.time.min(), ds.time.max()])
    ylim = ax1.get_ylim()
    ax1.set_ylim([ylim[0] * 0.5, ylim[1]])

    # title, labels, legend
    ax1.set_title('Flowline model')
    ax1.set_xlabel('Years of model evolution')
    leg_kwargs = {'bbox_to_anchor': (0.5, 0), 'loc': 8, 'ncol': 3}
    ax1.legend(handles_fl, labels_fl, title='Temperature bias', **leg_kwargs)
    ax1.yaxis.tick_right()
    ax1.yaxis.set_label_position('right')
    if ds.normalized:
        # add ylabel
        ax1.set_ylabel('Relative {}'.format(var))
        # aux line
        ax1.axhline(1, lw=0.8, ls=':', c='k')
    else:
        # add ylabel
        unit = 'm' if var == 'length' else 'm$^2$' if var == 'area' else 'm$^3$'
        ax1.set_ylabel('Glacier {} [{}]'.format(var, unit))
        # aux line
        ax1.axhline(ds.sel(model='fl')[var].isel(time=0).mean(), lw=0.8,
                    ls=':', c='k')

    ax1.grid()

    # add suptitle
    if suptitle:
        fig.suptitle(suptitle, fontsize=15)


def plot_both_climates(ds, var,
                       title=['Volume/area scaling model', 'Flowline model'],
                       xlim=None, scale=1):
    # define fig size and axes size and location
    figsize = [4, 3]
    ax_size = [0, 0, 1, 1]

    # ----------------------------
    #        SCALING MODEL
    # ----------------------------
    fig_vas = plt.figure(figsize=figsize)
    ax_vas = fig_vas.add_axes(ax_size)
    # flowline model
    (ds.sel(model='fl', mb_model='random')[var] * scale).plot(hue='temp_bias',
                                                              ax=ax_vas,
                                                              add_legend=False,
                                                              color='lightgray',
                                                              lw=0.5)
    # vas model
    ax_vas.set_prop_cycle('color', vas_cycle)
    handles_vas_r = (ds.sel(model='vas', mb_model='random')[var] * scale).plot(
        hue='temp_bias', ax=ax_vas, add_legend=False, lw=2)
    handles_vas_c = (
            ds.sel(model='vas', mb_model='constant')[var] * scale).plot(
        hue='temp_bias', ax=ax_vas, add_legend=False,
        lw=1.7, ls='--')
    # define legend labels
    labels_vas = ["$\\bf{const:}$", "$\\bf{rand:}$"]
    labels_vas.extend(np.array(
        ['{:+.1f} °C'.format(bias) for bias in ds.temp_bias.values]).repeat(2))
    # define legend handles
    title_proxy, = plt.plot(
        (ds.sel(model='vas', mb_model='random')[var] * scale).isel(temp_bias=0,
                                                                   time=0).values,
        marker='None', linestyle='None', label='dummy')
    handles_vas = list(np.array(np.repeat(title_proxy, 2)))
    handles_vas.extend(np.array([handles_vas_c, handles_vas_r]).T.flatten())

    # set axes limits
    if xlim:
        ax_vas.set_xlim(xlim)
    else:
        ax_vas.set_xlim([ds.time.min(), ds.time.max()])

    # title, labels, legend
    if title:
        ax_vas.set_title(title[0])
    else:
        ax_vas.set_title('')
    ax_vas.set_xlabel('Years of model evolution')
    ax_vas.legend(handles_vas, labels_vas, title='', bbox_to_anchor=(0.5, 1),
                  loc=8, ncol=4)
    if ds.normalized:
        # add ylabel
        ax_vas.set_ylabel('Relative {}'.format(var))
        # aux line
        ax_vas.axhline(1, lw=0.8, ls=':', c='k')
    else:
        # add ylabel
        unit = 'km' if var == 'length' else 'km$^2$' if var == 'area' else 'km$^3$'
        ax_vas.set_ylabel('Glacier {} [{}]'.format(var, unit))
        # aux line
        ax_vas.axhline((ds.sel(model='vas')[var] * scale).isel(time=0).mean(),
                       lw=0.8,
                       ls=':', c='k')

    # add grid
    ax_vas.grid()

    # ----------------------------
    #        FLOWLINE MODEL
    # ----------------------------
    fig_fl = plt.figure(figsize=figsize)
    ax_fl = fig_fl.add_axes(ax_size)

    # vas model
    (ds.sel(model='vas', mb_model='random')[var] * scale).plot(hue='temp_bias',
                                                               ax=ax_fl,
                                                               add_legend=False,
                                                               color='lightgray',
                                                               lw=0.5)
    # flowline model
    ax_fl.set_prop_cycle('color', fl_cycle)
    handles_fl_r = (ds.sel(model='fl', mb_model='random')[var] * scale).plot(
        hue='temp_bias', ax=ax_fl, add_legend=False, lw=2)
    handles_fl_c = (ds.sel(model='fl', mb_model='constant')[var] * scale).plot(
        hue='temp_bias', ax=ax_fl, add_legend=False,
        lw=1.7, ls='--')

    # define legend labels
    labels_fl = ["$\\bf{const:}$", "$\\bf{rand:}$"]
    labels_fl.extend(np.array(
        ['{:+.1f} °C'.format(bias) for bias in ds.temp_bias.values]).repeat(2))
    # define legend handles
    handles_fl = list(np.array(np.repeat(title_proxy, 2)))
    handles_fl.extend(np.array([handles_fl_c, handles_fl_r]).T.flatten())

    # set axes limits
    if xlim:
        ax_fl.set_xlim(xlim)
    else:
        ax_fl.set_xlim([ds.time.min(), ds.time.max()])

    # title, labels, legend
    if title:
        ax_fl.set_title(title[0])
    else:
        ax_fl.set_title('')
    ax_fl.set_xlabel('Years of model evolution')
    ax_fl.legend(handles_fl, labels_fl, bbox_to_anchor=(0.5, 1.0), loc=8,
                 ncol=4)
    ax_fl.yaxis.tick_right()
    ax_fl.yaxis.set_label_position('right')
    if ds.normalized:
        # add ylabel
        ax_fl.set_ylabel('Relative {}'.format(var))
        # aux line
        ax_fl.axhline(1, lw=0.8, ls=':', c='k')
    else:
        # add ylabel
        unit = 'km' if var == 'length' else 'km$^2$' if var == 'area' else 'km$^3$'
        ax_fl.set_ylabel('Glacier {} [{}]'.format(var, unit))
        # aux line
        ax_fl.axhline(ds.sel(model='fl')[var].isel(time=0).mean() * scale,
                      lw=0.8,
                      ls=':', c='k')

    # add grid
    ax_fl.grid()

    # assure same y-scale for both figures
    ylim_vas = ax_vas.get_ylim()
    ylim_fl = ax_fl.get_ylim()
    ylim = [np.min([ylim_vas, ylim_fl]), np.max([ylim_vas, ylim_fl])]
    ax_vas.set_ylim(ylim)
    ax_fl.set_ylim(ylim)

    return fig_vas, fig_fl


def plot_one_fig_per_model(ds, var,
                           title=['Volume/area scaling model',
                                  'Flowline model'],
                           xlim=None):
    # define fig size and axes size and location
    figsize = [4, 3]

    # ----------------------------
    #        SCALING MODEL
    # ----------------------------
    fig_vas = plt.figure(figsize=figsize)
    ax_size = [0.2, 0.2, 0.7, 0.7]
    ax_vas = fig_vas.add_axes(ax_size)
    # flowline model
    ds.sel(model='fl', mb_model='constant')[var].plot(hue='temp_bias',
                                                      ax=ax_vas,
                                                      add_legend=False,
                                                      color='lightgray',
                                                      lw=0.5)
    # vas model
    ax_vas.set_prop_cycle('color', vas_cycle)
    handles_vas = ds.sel(model='vas', mb_model='constant')[var].plot(
        hue='temp_bias', ax=ax_vas, add_legend=False, lw=2)
    # define legend labels
    labels_vas = np.array(
        ['{:+.1f} °C'.format(bias) for bias in ds.temp_bias.values])

    # set axes limits
    if xlim:
        ax_vas.set_xlim(xlim)
    else:
        ax_vas.set_xlim([ds.time.min(), ds.time.max()])

    # title, labels, legend
    if title:
        ax_vas.set_title(title[0])
    else:
        ax_vas.set_title('')
    ax_vas.set_xlabel('Years of model evolution')
    legend_kwargs = {'title': '', 'bbox_to_anchor': (0.5, 0.98),
                     'loc': 'upper center', 'ncol': 2}
    ax_vas.legend(handles_vas, labels_vas, **legend_kwargs)
    if ds.normalized:
        # add ylabel
        ax_vas.set_ylabel('Relative {}'.format(var))
        # aux line
        ax_vas.axhline(1, lw=0.8, ls=':', c='k')
    else:
        # add ylabel
        unit = 'm' if var == 'length' else 'm$^2$' if var == 'area' else 'm$^3$'
        ax_vas.set_ylabel('Glacier {} [{}]'.format(var, unit))
        # aux line
        ax_vas.axhline(ds.sel(model='vas')[var].isel(time=0).mean(), lw=0.8,
                       ls=':', c='k')

    # add grid
    ax_vas.grid()

    # ----------------------------
    #        FLOWLINE MODEL
    # ----------------------------
    fig_fl = plt.figure(figsize=figsize)
    ax_size = [0.1, 0.2, 0.7, 0.7]
    ax_fl = fig_fl.add_axes(ax_size)

    # vas model
    ds.sel(model='vas', mb_model='constant')[var].plot(hue='temp_bias',
                                                       ax=ax_fl,
                                                       add_legend=False,
                                                       color='lightgray',
                                                       lw=0.5)
    # flowline model
    ax_fl.set_prop_cycle('color', fl_cycle)
    handles_fl = ds.sel(model='fl', mb_model='constant')[var].plot(
        hue='temp_bias', ax=ax_fl, add_legend=False, lw=2)

    # define legend labels
    labels_fl = np.array(
        ['{:+.1f} °C'.format(bias) for bias in ds.temp_bias.values])

    # set axes limits
    if xlim:
        ax_fl.set_xlim(xlim)
    else:
        ax_fl.set_xlim([ds.time.min(), ds.time.max()])

    # title, labels, legend
    if title:
        ax_fl.set_title(title[0])
    else:
        ax_fl.set_title('')
    ax_fl.set_xlabel('Years of model evolution')
    ax_fl.legend(handles_fl, labels_fl, **legend_kwargs)
    ax_fl.yaxis.tick_right()
    ax_fl.yaxis.set_label_position('right')
    if ds.normalized:
        # add ylabel
        ax_fl.set_ylabel('Relative {}'.format(var))
        # aux line
        ax_fl.axhline(1, lw=0.8, ls=':', c='k')
    else:
        # add ylabel
        unit = 'm' if var == 'length' else 'm$^2$' if var == 'area' else 'm$^3$'
        ax_fl.set_ylabel('Glacier {} [{}]'.format(var, unit))
        # aux line
        ax_fl.axhline(ds.sel(model='fl')[var].isel(time=0).mean(), lw=0.8,
                      ls=':', c='k')

    # add grid
    ax_fl.grid()

    # assure same y-scale for both figures
    ylim_vas = ax_vas.get_ylim()
    ylim_fl = ax_fl.get_ylim()
    ylim = [np.min([ylim_vas, ylim_fl]), np.max([ylim_vas, ylim_fl])]
    ax_vas.set_ylim(ylim)
    ax_fl.set_ylim(ylim)

    return fig_vas, fig_fl


def plot_both_climates_same_fig(ds, var, title='', suptitle='', xlim=None):
    """

    Parameters
    ----------
    ds
    var
    title
    suptitle
    xlim

    """
    # plot relative volume change
    fig, [ax0, ax1] = plt.subplots(1, 2, figsize=[12, 5])

    # SCALING MODEL - right panel
    # ---------------------------

    # flowline model
    ds.sel(model='fl', mb_model='constant')[var].plot(hue='temp_bias', ax=ax0,
                                                      add_legend=False,
                                                      color='lightgray',
                                                      lw=0.5)
    # vas model
    ax0.set_prop_cycle('color', vas_cycle)
    handles_vas_r = ds.sel(model='vas', mb_model='random')[var].plot(
        hue='temp_bias', ax=ax0, add_legend=False, lw=2)
    handles_vas_c = ds.sel(model='vas', mb_model='constant')[var].plot(
        hue='temp_bias', ax=ax0, add_legend=False,
        lw=1.7, ls='--')

    # define legend labels
    labels_vas = ["$\\bf{const:}$", "$\\bf{rand:}$"]
    labels_vas.extend(np.array(
        ['{:+.1f} °C'.format(bias) for bias in ds.temp_bias.values]).repeat(2))
    # define legend handles
    title_proxy, = plt.plot(0, marker='None', linestyle='None', label='dummy')
    handles_vas = list(np.array(np.repeat(title_proxy, 2)))
    handles_vas.extend(np.array([handles_vas_c, handles_vas_r]).T.flatten())

    # set axes limits
    if xlim:
        ax0.set_xlim(xlim)
    else:
        ax0.set_xlim([ds.time.min(), ds.time.max()])
    ylim = ax0.get_ylim()
    ax0.set_ylim([ylim[0] * 0.5, ylim[1]])

    # title, labels, legend
    ax0.set_title('Volume/area scaling model')
    ax0.set_xlabel('Years of model evolution')
    ax0.legend(handles_vas, labels_vas, title='',
               bbox_to_anchor=(0.5, 0), loc=8, ncol=4)
    if ds.normalized:
        # add ylabel
        ax0.set_ylabel('Relative {}'.format(var))
        # aux line
        ax0.axhline(1, lw=0.8, ls=':', c='k')
    else:
        # add ylabel
        unit = 'm' if var == 'length' else 'm$^2$' if var == 'area' else 'm$^3$'
        ax0.set_ylabel('Glacier {} [{}]'.format(var, unit))
        # aux line
        ax0.axhline(ds.sel(model='vas')[var].isel(time=0).mean(), lw=0.8,
                    ls=':', c='k')

    ax0.grid()

    # FLOWLINE MODEL - right panel
    # ----------------------------

    # vas model
    ds.sel(model='vas', mb_model='constant')[var].plot(hue='temp_bias', ax=ax1,
                                                       add_legend=False,
                                                       color='lightgray',
                                                       lw=0.5)
    # flowline model
    ax1.set_prop_cycle('color', fl_cycle)
    handles_fl_r = ds.sel(model='fl', mb_model='random')[var].plot(
        hue='temp_bias', ax=ax1, add_legend=False, lw=2)
    handles_fl_c = ds.sel(model='fl', mb_model='constant')[var].plot(
        hue='temp_bias', ax=ax1, add_legend=False, lw=1.7,
        ls='--')

    # define legend labels
    labels_fl = ["$\\bf{const:}$", "$\\bf{rand:}$"]
    labels_fl.extend(np.array(
        ['{:+.1f} °C'.format(bias) for bias in ds.temp_bias.values]).repeat(2))
    # define legend handles
    handles_fl = list(np.array(np.repeat(title_proxy, 2)))
    handles_fl.extend(np.array([handles_fl_c, handles_fl_r]).T.flatten())

    # set axes limits
    if xlim:
        ax1.set_xlim(xlim)
    else:
        ax1.set_xlim([ds.time.min(), ds.time.max()])
    ylim = ax1.get_ylim()
    ax1.set_ylim([ylim[0] * 0.5, ylim[1]])

    # title, labels, legend
    ax1.set_title('Flowline model')
    ax1.set_xlabel('Years of model evolution')
    ax1.legend(handles_fl, labels_fl, bbox_to_anchor=(0.5, 0), loc=8, ncol=4)
    ax1.yaxis.tick_right()
    ax1.yaxis.set_label_position('right')
    if ds.normalized:
        # add ylabel
        ax1.set_ylabel('Relative {}'.format(var))
        # aux line
        ax1.axhline(1, lw=0.8, ls=':', c='k')
    else:
        # add ylabel
        unit = 'm' if var == 'length' else 'm$^2$' if var == 'area' else 'm$^3$'
        ax1.set_ylabel('Glacier {} [{}]'.format(var, unit))
        # aux line
        ax1.axhline(ds.sel(model='fl')[var].isel(time=0).mean(), lw=0.8,
                    ls=':', c='k')

    ax1.grid()

    # add suptitle
    if suptitle:
        fig.suptitle(suptitle, fontsize=15)


def plot_deviations_three_panel(ds, rgi_df, var='length', time0=None,
                                time1=None, out_dir=''):
    # iterate over all glaciers
    for i, (rgi_id, glacier) in enumerate(rgi_df.iterrows()):

        # create new figure and axes for each glacier
        fig, axes = plt.subplots(3, 1)

        # create empty containers
        ylim = list()
        ylim_ = list()
        axes_ = list()
        handles_fl = list()
        handles_vas = list()

        # iterate over all temperature biases
        for j, (ax, temp_bias) in enumerate(zip(axes, ds.temp_bias)):
            # select glacier and climate
            ds_var = ds[var].sel(rgi_id=rgi_id,
                                 temp_bias=temp_bias,
                                 normalized=False).squeeze()
            # select time period
            time0 = (int(8e3) if not time0 else time0)
            time1 = (int(12e3) + 1 if not time1 else time1)
            ds_var = ds_var.isel(time=slice(time0, time1))

            # compute mean and deviations
            mean = ds_var.mean(dim='time')
            deviations = ds_var - mean

            # plot var anomalies
            handles_vas.append(ax.plot(deviations.time,
                                       deviations.sel(model='vas'),
                                       c=vas_cycle[j], lw=2,
                                       label=f'vas {temp_bias.values:+.1f}'))

            # instantiate a second axes that shares the same x-axis
            ax_ = ax.twinx()
            axes_.append(ax_)

            handles_fl.append(ax_.plot(deviations.time,
                                       deviations.sel(model='fl'),
                                       c=fl_cycle[j],
                                       label=f'fl {temp_bias.values:+.1f}'))
            ax_.axhline(0, c='k', ls=':', lw=0.8)

            ylim.append(ax.get_ylim())
            ylim_.append(ax_.get_ylim())

            # take care of tickparams
            ax.tick_params(direction='inout', top=True, labelbottom=(j == 2))
            # set xlimits to match sliced time
            ax.set_xlim([time0, time1])

        # title, labels
        fig.text(0.5, 0.02, 'Years of model evolution', ha='center')
        unit = 'm$^3$' if var == 'volume' else (
            'm$^2$' if var == 'area' else 'm')
        fig.text(1.02, 0.5, f'Flowline {var} anomalies [{unit}]', va='center',
                 rotation='vertical')
        fig.text(-0.02, 0.5, f'VAS {var} anomalies [{unit}]', va='center',
                 rotation='vertical')

        # adjust layout
        fig.subplots_adjust(wspace=0, hspace=0)

        # adjust limits
        ylim = abs(np.array(ylim).flatten()).max()
        [ax.set_ylim([-ylim, ylim]) for ax in axes]
        ylim_ = abs(np.array(ylim_).flatten()).max()
        [ax_.set_ylim([-ylim_, ylim_]) for ax_ in axes_]

        # define legend handles and labels
        title_proxy, = plt.plot(0, time0, marker='None', linestyle='None',
                                label='dummy')
        # define legend labels
        labels = ["$\\bf{VAS:}$", "$\\bf{Flowline:}$"]
        labels.extend(np.array(
            ['{:+.1f} °C'.format(bias) for bias in
             ds.temp_bias.values]).repeat(2))
        handles = list(np.array(np.repeat(title_proxy, 2)))
        handles.extend(
            np.array([handles_vas, handles_fl]).T.flatten())
        fig.legend(handles, labels, bbox_to_anchor=(0.5, 1.0),
                   loc='upper center',
                   ncol=4)

        # save to file
        f_name = f"{glacier['name'].replace(' ', '_')}.pdf"
        fig.savefig(os.path.join(out_dir, f_name), bbox_inches='tight')


def plot_deviations_three_panel_old(ds, rgi_df, var='length', time0=None,
                                    time1=None, out_dir=''):
    # iterate over all glaciers
    for i, (rgi_id, glacier) in enumerate(rgi_df.iterrows()):

        # create new figure and axes for each glacier
        fig, axes = plt.subplots(3, 1)

        # create empty containers
        ylim = list()
        ylim_ = list()
        axes_ = list()
        handles_fl = list()
        handles_vas = list()

        # iterate over all temperature biases
        for j, (ax, temp_bias) in enumerate(zip(axes, ds.temp_bias)):
            # select glacier and climate
            ds_var = ds[var].sel(rgi_id=rgi_id,
                                 temp_bias=temp_bias,
                                 normalized=False).squeeze()
            # select time period
            time0 = (int(8e3) if not time0 else time0)
            time1 = (int(12e3) + 1 if not time1 else time1)
            ds_var = ds_var.isel(time=slice(time0, time1))

            # compute mean and deviations
            mean = ds_var.mean(dim='time')
            deviations = ds_var - mean

            # plot var anomalies
            handles_fl.append(ax.plot(deviations.time,
                                      deviations.sel(model='fl'),
                                      c=fl_cycle[j],
                                      label=f'fl {temp_bias.values:+.1f}'))

            # instantiate a second axes that shares the same x-axis
            ax_ = ax.twinx()
            axes_.append(ax_)
            handles_vas.append(ax_.plot(deviations.time,
                                        deviations.sel(model='vas'),
                                        c=vas_cycle[j], lw=2,
                                        label=f'vas {temp_bias.values:+.1f}'))
            ax_.axhline(0, c='k', ls=':', lw=0.8)

            ylim.append(ax.get_ylim())
            ylim_.append(ax_.get_ylim())

            # take care of tickparams
            ax.tick_params(direction='inout', top=True, labelbottom=(j == 2))
            # set xlimits to match sliced time
            ax.set_xlim([time0, time1])

        # title, labels
        fig.text(0.5, 0.02, 'Years of model evolution', ha='center')
        unit = 'm$^3$' if var == 'volume' else (
            'm$^2$' if var == 'area' else 'm')
        fig.text(-0.02, 0.5, f'Flowline {var} anomalies [{unit}]', va='center',
                 rotation='vertical')
        fig.text(1.02, 0.5, f'VAS {var} anomalies [{unit}]', va='center',
                 rotation='vertical')

        # adjust layout
        fig.subplots_adjust(wspace=0, hspace=0)

        # adjust limits
        ylim = abs(np.array(ylim).flatten()).max()
        [ax.set_ylim([-ylim, ylim]) for ax in axes]
        ylim_ = abs(np.array(ylim_).flatten()).max()
        [ax_.set_ylim([-ylim_, ylim_]) for ax_ in axes_]

        # define legend handles and labels
        title_proxy, = plt.plot(0, time0, marker='None', linestyle='None',
                                label='dummy')
        # define legend labels
        labels = ["$\\bf{VAS:}$", "$\\bf{Flowline:}$"]
        labels.extend(np.array(
            ['{:+.1f} °C'.format(bias) for bias in
             ds.temp_bias.values]).repeat(2))
        handles = list(np.array(np.repeat(title_proxy, 2)))
        handles.extend(
            np.array([handles_vas, handles_fl]).T.flatten())
        fig.legend(handles, labels, bbox_to_anchor=(0.5, 1.0),
                   loc='upper center',
                   ncol=4)

        # save to file
        f_name = f"{glacier['name'].replace(' ', '_')}.pdf"
        fig.savefig(os.path.join(out_dir, f_name), bbox_inches='tight')


# --------------------------------------------
# Calling the plotting functions on real data
# --------------------------------------------


def plot_showcase_glaciers_random_climate():
    # specify path and read datasets
    path = '/Users/oberrauch/work/master/data/cluster_output/showcase_glaciers_random_climate/eq_runs.nc'
    ds = read_dataset(path)

    # define glaciers of interest
    showcase_glaciers = pd.read_csv(
        '/Users/oberrauch/work/master/data/showcase_glaciers.csv', index_col=0)

    # iterate over all above selected glaciers
    for rgi_id, glacier in showcase_glaciers.iterrows():
        # select glacier
        rgi_id = rgi_id
        name = glacier['name']
        log.info('Plots for {} ({})'.format(name, rgi_id))

        # iterate over all selected variables
        variables = ['length', 'volume']
        for var in variables:
            for norm in [True, False]:
                log.info('Plotting {} {} time series'
                         .format('normalized' if norm else 'absolute',
                                 var))
                suptitle = '{} evolution under random climate'.format(name)
                plot_eq_time_series(ds.sel(mb_model='random',
                                           normalized=norm,
                                           rgi_id=rgi_id),
                                    var=var, suptitle=suptitle, xlim=[0, 1e3])
                f_path = '/Users/oberrauch/work/master/plots/final_plots/time_series/showcase_glaciers_random_climate/'
                f_path += '{}_{}_{}.pdf'.format(var, 'norm' if norm else 'abs',
                                                name.replace(' ', '_'))
                plt.savefig(f_path, bbox_inches='tight')
                plt.close(plt.gcf())

    # close dataset
    ds.close()


def plot_single_glaciers_both_climates():
    # specify path and read datasets
    path = '/Users/oberrauch/work/master/data/cluster_output/single_glaciers/eq_runs.nc'
    ds = read_dataset(path)

    # define glaciers of interest
    showcase_glaciers = pd.read_csv(
        '/Users/oberrauch/work/master/data/showcase_glaciers.csv', index_col=0)

    # iterate over all above selected glaciers
    for rgi_id, glacier in showcase_glaciers.iterrows():
        # select glacier
        rgi_id = rgi_id
        name = glacier['name']
        log.info('Plots for {} ({})'.format(name, rgi_id))

        # iterate over all selected variables
        variables = ['length', 'volume', 'area']
        for var in variables:
            # iterate over normalized and absolute values
            # for norm in [True, False]:
            for norm in [True, False]:
                log.info('Plotting {} {} time series'
                         .format('normalized' if norm else 'absolute',
                                 var))
                # define limit of x-axis
                xlim = [0, 1e3]

                figs = plot_both_climates(
                    ds.sel(normalized=norm, rgi_id=rgi_id),
                    var=var, title=None, xlim=xlim)

                # store to file
                for m, fig in zip(['vas', 'fl'], figs):
                    f_path = '/Users/oberrauch/work/master/plots/final_plots/time_series/single_glaciers/'
                    f_path += '{}_{}_{}_{}.pdf'.format(var,
                                                       'norm' if norm else 'abs',
                                                       m,
                                                       name.replace(' ', '_'))
                    fig.savefig(f_path, bbox_inches='tight')
                    plt.close(fig)

    # close dataset
    ds.close()


def plot_histalp_commitment():
    # specify path and read datasets
    path = '/Users/oberrauch/work/master/data/cluster_output/histalp_commitment/eq_runs.nc'
    ds = read_dataset(path)

    # iterate over mass balance model
    for mb_model in ds['mb_model'].values:
        log.info('Creating plots for {} climate scenario'.format(mb_model))
        # iterate over all selected variables
        variables = ['length', 'volume']
        for var in variables:
            # iterate over normalized and absolute values
            for norm in [True, False]:
                log.info('Plotting {} {} time series'
                         .format('normalized' if norm else 'absolute',
                                 var))
                # define limit of x-axis depending on climate scenario
                xlim = [0, 1e4] if mb_model == 'random' else [0, 1e3]

                suptitle = None
                plot_eq_time_series(
                    ds.sel(mb_model=mb_model, normalized=norm, rgi_id='sum'),
                    var=var, suptitle=suptitle, xlim=xlim)
                # store to file
                f_path = '/Users/oberrauch/work/master/plots/final_plots/time_series/histalp_commitment/'
                f_path += '{}_{}_{}.pdf'.format(var, 'norm' if norm else 'abs',
                                                mb_model)
                plt.savefig(f_path, bbox_inches='tight')
                plt.close(plt.gcf())

    # close dataset
    ds.close()


def plot_histalp_commitment_both_climates():
    """
    """
    # specify path and read dataset
    path = '/Users/oberrauch/work/master/data/cluster_output/histalp_commitment/eq_runs.nc'
    ds = read_dataset(path)

    # iterate over all selected variables
    variables = ['length', 'volume']
    for var in variables:
        # iterate over normalized and absolute values
        for norm in [True, False]:
            log.info('Plotting {} {} time series'
                     .format('normalized' if norm else 'absolute',
                             var))
            # define limit of x-axis
            xlim = [0, 1e3]

            scale = 1
            if not norm:
                if var == 'length':
                    scale = 1e-3
                if var == 'volume':
                    scale = 1e-9

            figs = plot_both_climates(ds.sel(normalized=norm, rgi_id='sum'),
                                      var=var, title=None, xlim=xlim,
                                      scale=scale)

            # store to file
            for m, fig in zip(['vas', 'fl'], figs):
                f_path = '/Users/oberrauch/work/master/plots/final_plots/time_series/histalp_commitment/'
                f_path += '{}_{}_{}.pdf'.format(var, 'norm' if norm else 'abs',
                                                m)
                fig.savefig(f_path, bbox_inches='tight')
                plt.close(fig)

    # close dataset
    ds.close()


def plot_histalp_projection():
    """
    """
    # specify path and read dataset
    path = '/Users/oberrauch/work/master/data/cluster_output/histalp_projection/eq_runs.nc'
    ds = read_dataset(path)
    ds = ds.sel(temp_bias=[0, 0.5, 1, 2])

    # iterate over all selected variables
    variables = ['length', 'volume']
    for var in variables:
        # iterate over normalized and absolute values
        for norm in [True, False]:
            log.info('Plotting {} {} time series'
                     .format('normalized' if norm else 'absolute',
                             var))
            # define limit of x-axis
            xlim = [0, 300]

            figs = plot_one_fig_per_model(
                ds.sel(normalized=norm, rgi_id='sum'),
                var=var, title=None, xlim=xlim)

            # store to file
            for m, fig in zip(['vas', 'fl'], figs):
                f_path = '/Users/oberrauch/work/master/plots/final_plots/time_series/histalp_projection/'
                f_path += '{}_{}_{}.pdf'.format(var, 'norm' if norm else 'abs',
                                                m)
                fig.savefig(f_path)
                plt.close(fig)

    # close dataset
    ds.close()


def plot_random_length():
    # specify path and read datasets
    path = '/Users/oberrauch/work/master/data/cluster_output/showcase_glaciers_random_climate/eq_runs.nc'
    ds = read_dataset(path)

    # define glaciers of interest
    showcase_glaciers = pd.read_csv(
        '/Users/oberrauch/work/master/data/showcase_glaciers.csv', index_col=0)

    out_dir = '/Users/oberrauch/work/master/plots/final_plots/random_length/'
    plot_deviations_three_panel(ds, rgi_df=showcase_glaciers, out_dir=out_dir)


if __name__ == '__main__':
    # plot_single_glaciers_both_climates()
    plot_random_length()
    # plot_histalp_commitment_both_climates()
    # plot_histalp_projection()
