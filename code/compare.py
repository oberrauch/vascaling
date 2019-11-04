""" Evolution of length, area and volume of HEF over the HistAlp climate
period (1802-2003). Comparing the volume/area scaling approach to the
OGGM flowline model dynamics.

Old: Data generated by `run_vas.py` and `run_oggm.py`, respectively.
Edit: use data generated by `run_comparison.py` instead (16/01/19).
"""

# import external libs
from os import path
import pandas as pd
import matplotlib.pyplot as plt

# load DataFrame
folder = '/Users/oberrauch/work/master/data/'
df = pd.read_csv(path.join(folder, 'run_comparison.csv'), index_col=0)


def plot_twin_x(vas, oggm, title='', ylabel='', file_path=''):
    """ Plot geometric parameters of both models,
        the result of the scaling model on the left y-axis and
        the result of the OGGM on the right y-axis.
        If a `file_path` is given, the figure will be saved.

    :param vas: (array like) geometric glacier parameter of the VAS model
    :param oggm: (array like) geometric glacier parameter of the flowline model
    :param title: (string) figure title
    :param ylabel: (string) label for y-axis
    :param file_path: (string) where to store the figure
    """
    # create figure and first axes
    fig, ax1 = plt.subplots()
    # define first color
    color = 'C0'
    # plot
    ax1.plot(vas, c=color)
    ax1.set_ylabel(ylabel + ' (Scaling)', color=color)
    ax1.tick_params(axis='y', labelcolor=color)
    # create second axes
    ax2 = ax1.twinx()
    # plot
    color = 'C1'
    ax2.plot(oggm, c=color)
    ax2.set_ylabel(ylabel + ' (OGGM)', color=color)
    ax2.tick_params(axis='y', labelcolor=color)
    # add title
    plt.title(title)

    # store to file
    if file_path:
        plt.savefig(file_path, bbox_inches='tight')


def plot_both(vas, oggm, ref=None, title='', ylabel='', file_path='', exp=0):
    """ Plot geometric parameters of both models.
    If a `file_path` is given, the figure will be saved.

    :param vas: (pandas.Series) geometric glacier parameter of the VAS model
    :param oggm: (pandas.Series) geometric glacier parameter of the OGGM
    :param ref: (pandas.Series) measured glacier parameter, optional
    :param title: (string) figure title, optional
    :param ylabel: (string) label for y-axis, optional
    :param file_path: (string) where to store the figure, optional
    :param exp: (int) exponent for labels in scientific notation, optional
    """
    # create figure and first axes
    fig = plt.figure(figsize=[6, 4])
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

    # plot
    c1 = 'C0'
    c2 = 'C1'
    c3 = 'C3'
    ax.plot(vas, c=c1, label='VAS')
    ax.plot(oggm, c=c2, label='OGGM')
    if ref:
        ax.plot(ref, c=c3, label='measurements')

    # add title, labels, legend
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.legend()
    # use scientific notation with fixed exponent according
    ax.ticklabel_format(style='sci', axis='y', scilimits=(exp, exp))

    # store to file
    if file_path and False:
        plt.savefig(file_path, bbox_inches='tight')


# specify plot directory
folder = '/Users/oberrauch/work/master/plots/'
glacier_name = 'Ob. Grindelwaldgletscher'

# plot length
plot_both(df.length_vas, df.length_oggm,
          title='{} Length'.format(glacier_name),
          ylabel='Length [km]',
          file_path=folder+'length.png',
          exp=3)
# plot area
plot_both(df.area_vas, df.area_oggm,
          title='{} area'.format(glacier_name),
          ylabel='Area [km2]',
          file_path=folder+'area.png',
          exp=6)
# plot volume
plot_both(df.volume_vas, df.volume_oggm,
          title='{} volume'.format(glacier_name),
          ylabel='Volume [km3]',
          file_path=folder+'volume.png',
          exp=9)
