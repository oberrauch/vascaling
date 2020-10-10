""" Color schemes for my masters thesis. Includes color cycles and different
getter-functions. Color maps generated with coolors.co.

"""
import numpy as np
import matplotlib.pyplot as plt

# Flowline color palette
# coolors.co/ddcb64-8db854-14381e-5a72ac-3d2c76-0b111e-edbad9-821753-411148

fl_cycle = np.array(["#57b9ff", "#1978b8", "#08102b", "#ddcb64", "#8db854", "#14381e", "#edbad9", "#821753", "#411148"])
fl_dict = {"Arylide Yellow": "#ddcb64", "Bud Green": "#8db854", "Phthalo Green": "#14381e", "Blue Yonder": "#5a72ac",
           "Spanish Violet": "#3d2c76", "Rich Black FOGRA 29": "#0b111e", "Cotton Candy": "#edbad9",
           "Pansy Purple": "#821753", "Russian Violet": "#411148"}
fl_extended = [
    {"name": "Arylide Yellow", "hex": "#ddcb64", "rgb": [221, 203, 100], "cmyk": [0, 8, 55, 13], "hsb": [51, 55, 87],
     "hsl": [51, 64, 63], "lab": [81, -7, 53]},
    {"name": "Bud Green", "hex": "#8db854", "rgb": [141, 184, 84], "cmyk": [23, 0, 54, 28], "hsb": [86, 54, 72],
     "hsl": [86, 41, 53], "lab": [70, -31, 46]},
    {"name": "Phthalo Green", "hex": "#14381e", "rgb": [20, 56, 30], "cmyk": [64, 0, 46, 78], "hsb": [137, 64, 22],
     "hsl": [137, 47, 15], "lab": [20, -20, 12]},
    {"name": "Blue Yonder", "hex": "#5a72ac", "rgb": [90, 114, 172], "cmyk": [48, 34, 0, 33], "hsb": [222, 48, 67],
     "hsl": [222, 33, 51], "lab": [48, 7, -34]},
    {"name": "Spanish Violet", "hex": "#3d2c76", "rgb": [61, 44, 118], "cmyk": [48, 63, 0, 54], "hsb": [254, 63, 46],
     "hsl": [254, 46, 32], "lab": [24, 28, -40]},
    {"name": "Rich Black FOGRA 29", "hex": "#0b111e", "rgb": [11, 17, 30], "cmyk": [63, 43, 0, 88],
     "hsb": [221, 63, 12],
     "hsl": [221, 46, 8], "lab": [5, 1, -9]},
    {"name": "Cotton Candy", "hex": "#edbad9", "rgb": [237, 186, 217], "cmyk": [0, 22, 8, 7], "hsb": [324, 22, 93],
     "hsl": [324, 59, 83], "lab": [81, 23, -8]},
    {"name": "Pansy Purple", "hex": "#821753", "rgb": [130, 23, 83], "cmyk": [0, 82, 36, 49], "hsb": [326, 82, 51],
     "hsl": [326, 70, 30], "lab": [29, 49, -8]},
    {"name": "Russian Violet", "hex": "#411148", "rgb": [65, 17, 72], "cmyk": [10, 76, 0, 72], "hsb": [292, 76, 28],
     "hsl": [292, 62, 17], "lab": [15, 32, -23]}]

# Volume/area scaling palette
# coolors.co/f9cb16-f69d0e-ca3221-57b9ff-1978b8-08102b-ffa8aa-b788c8-31276d

vas_cycle = np.array(
    ["#f9cb16", "#f69d0e", "#ca3221", "#5a72ac", "#3d2c76", "#0b111e", "#ffa8aa", "#b788c8", "#31276d"])

vas_dict = {"Jonquil": "#f9cb16", "Orange Peel": "#f69d0e", "Maximum Red": "#ca3221", "Maya Blue": "#57b9ff",
            "Star Command Blue": "#1978b8", "Oxford Blue": "#08102b", "Light Pink": "#ffa8aa",
            "African Violet": "#b788c8",
            "St Patricks Blue": "#31276d"}
vas_extended = [
    {"name": "Jonquil", "hex": "#f9cb16", "rgb": [249, 203, 22], "cmyk": [0, 18, 91, 2], "hsb": [48, 91, 98],
     "hsl": [48, 95, 53], "lab": [83, 2, 82]},
    {"name": "Orange Peel", "hex": "#f69d0e", "rgb": [246, 157, 14], "cmyk": [0, 36, 94, 4],
     "hsb": [37, 94, 96], "hsl": [37, 93, 51], "lab": [72, 24, 75]},
    {"name": "Maximum Red", "hex": "#ca3221", "rgb": [202, 50, 33], "cmyk": [0, 75, 84, 21],
     "hsb": [6, 84, 79], "hsl": [6, 72, 46], "lab": [46, 58, 46]},
    {"name": "Maya Blue", "hex": "#57b9ff", "rgb": [87, 185, 255], "cmyk": [66, 27, 0, 0],
     "hsb": [205, 66, 100], "hsl": [205, 100, 67], "lab": [72, -8, -43]},
    {"name": "Star Command Blue", "hex": "#1978b8", "rgb": [25, 120, 184], "cmyk": [86, 35, 0, 28],
     "hsb": [204, 86, 72], "hsl": [204, 76, 41], "lab": [48, -3, -41]},
    {"name": "Oxford Blue", "hex": "#08102b", "rgb": [8, 16, 43], "cmyk": [81, 63, 0, 83],
     "hsb": [226, 81, 17], "hsl": [226, 69, 10], "lab": [5, 6, -19]},
    {"name": "Light Pink", "hex": "#ffa8aa", "rgb": [255, 168, 170], "cmyk": [0, 34, 33, 0],
     "hsb": [359, 34, 100], "hsl": [359, 100, 83], "lab": [77, 32, 12]},
    {"name": "African Violet", "hex": "#b788c8", "rgb": [183, 136, 200], "cmyk": [8, 32, 0, 22],
     "hsb": [284, 32, 78], "hsl": [284, 37, 66], "lab": [63, 29, -26]},
    {"name": "St Patricks Blue", "hex": "#31276d", "rgb": [49, 39, 109], "cmyk": [55, 64, 0, 57],
     "hsb": [249, 64, 43], "hsl": [249, 47, 29], "lab": [21, 25, -39]}]


def plot_colors(color_cylces, titles=None):
    """

    """
    # create figure and axes
    n_cycles = len(color_cylces)
    fig, ax = plt.subplots(1, n_cycles)

    # create empty titles if none are given
    if not titles:
        titles = np.repeat('', n_cycles)

    # iterate over all given color cycles
    for i, (color_cycle, title) in enumerate(zip(color_cylces, titles)):
        # iterate over all colos
        for j, c in enumerate(color_cycle):
            # plot horizontal line
            ax[i].axhline(j, c=c, lw='25')

        # add title
        ax[i].set_title(title, fontsize=18)

    # turn off axis
    for ax_ in ax:
        ax_.invert_yaxis()
        ax_.axis('off')

    plt.show()


if __name__ == '__main__':
    plot_colors([vas_cycle, fl_cycle], ["VAS", "Flowline"])
    pass
