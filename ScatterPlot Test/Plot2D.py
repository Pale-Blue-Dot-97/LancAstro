# =====================================================================================================================
# ===================================================== PLOT2D ========================================================
# =====================================================================================================================
"""Script for 2D plotting of data stored in FITS files, ascii tables

Part of the LancAstro.py project for XGAL

Author: Harry Baker
Date: 31.07.2019
Version: 1.3
== RELEASE ==

Example:
    $ SCPTest2.py

Todo:
    * Move parameters out to a separate config file

"""

# =====================================================================================================================
#                                                      IMPORTS
# =====================================================================================================================

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import random
import Dependencies.progress_bar as pb


# =====================================================================================================================
#                                                     FUNCTIONS
# =====================================================================================================================


def convert_to_array(var):
    """Checks if parameter supplied is an array and corrects if not

    Args:
        var (np.array/array/str/float): Parameter value/array to be checked and corrected

    Return:
        var
    """

    if not isinstance(var, (list, tuple, np.ndarray)):
        var = np.array(var)

    return var


def default_lengths(param, n_var):
    """
    Checks every parameter list is the same length. i.e if user only defines min and max for 1 variable.
    Corrects these to all have same length 

    Args:
        param (np.array): Parameter array to be set to correct length
        n_var (int): Number of variables therefore the length to set param to

    Return:
        param 

    """

    # Uses convert_to_array to make sure entered params are arrays
    param = convert_to_array(param)

    # Checks if every parameter is the length of the number of variables
    if len(param) != n_var:
        if print_on:
            print("WARNING: Not every variable has defined range limits!")
            print("Setting all variables to default of 1st entry \n")

        for i in range(len(param)):
            param[i] = param[0]

        for i in range(n_var - len(param)):
            param.append(param[0])

    return param


def colour_cycle(n_var, COLOURS):
    """Checks if a colour is defined for every variable and if not, appends the missing entries to random RGBa values

    Args:
        n_var (int): Number of variables to be plotted
        COLOURS ([str]): List of strings defining colours for each plot

    Returns:
        COLOURS ([str]): Amended list of colour definitions

    """

    if len(COLOURS) != n_var:
        for i in range(n_var):
            r = lambda: random.randint(0, 255)
            random_colour = '#%02X%02X%02X' % (r(), r(), r())
            if i < len(COLOURS):
                COLOURS[i] = random_colour
            else:
                COLOURS.append(random_colour)


def construct_ticks(axis_range, num_x, num_y):
    """Creates lists defining the positions of x and y axis ticks

    Args:
        axis_range (np.array): Range of axis
        num_x (int): Number of x-ticks desired
        num_y (int): Number of y-ticks desired

    Returns:
        x_ticks ([str]): Positions of x-axis ticks
        y_ticks ([str]): Positions of y-axis ticks

    """

    x_step = (axis_range[1] - axis_range[0]) / float(num_x)
    y_step = (axis_range[3] - axis_range[2]) / float(num_y)

    x_ticks = np.round(np.arange(axis_range[0], axis_range[1], x_step), 1)
    y_ticks = np.round(np.arange(axis_range[2], axis_range[3], y_step), 1)

    return x_ticks, y_ticks


def plot(x, y, x_error=[None], y_error=[None], DATALABELS=[None], COLOURS=[None], FILL_COLOURS=[None], SHADE=[0.5],
         SIZES=[None], POINTSTYLES=['none'], EDGECOLOURS=[None], EDGEWID=[None], LWID=[2], ERBWID=[2]):
    """ Method capable of plotting multiple sets of 2D variables

    Args:
        x ([[float]]): All arrays of x-axis data to plot
        y ([[float]]): All arrays of y-axis data to plot
        x_error ([[float]]/[[tuple]]): Errors for each x value. Supports tuples for asymmetric errors
        y_error ([[float]]/[[tuple]]): Errors for each y value. Supports tuples for asymmetric errors
        DATALABELS ([str]): Labels for legend of plots
        COLOURS ([str]): Colours for each plot
        FILL_COLOURS ([str]): Colour of fill for plot. If None is entered then no fill is used
        SHADE ([float]): Alpha value for fill under each plot
        SIZES ([float]): Size of points for each plot
        POINTSTYLES ([str]): Style of point/ line for each plot
        EDGECOLOURS ([str]): Colour for points of each plot
        EDGEWID ([float]): Width of point edge for each plot
        LWID ([float]): Width of line for each plot
        ERBWID ([float]): Width of error bars for each plot

    Returns:
        HANDLES: Handles for each plot to be used to construct a legend and final figure
    """

    HANDLES = []

    x = convert_to_array(x)
    y = convert_to_array(y)

    n_var = len(x)

    # Makes sure all parameters are the same length as number of plots
    x_error = default_lengths(x_error, n_var)
    y_error = default_lengths(y_error, n_var)
    DATALABELS = default_lengths(DATALABELS, n_var)
    COLOURS = default_lengths(COLOURS, n_var)
    FILL_COLOURS = default_lengths(FILL_COLOURS, n_var)
    SIZES = default_lengths(SIZES, n_var)
    POINTSTYLES = default_lengths(POINTSTYLES, n_var)
    EDGECOLOURS = default_lengths(EDGECOLOURS, n_var)
    SHADE = default_lengths(SHADE, n_var)
    LWID = default_lengths(LWID, n_var)
    ERBWID = default_lengths(ERBWID, n_var)
    EDGEWID = default_lengths(EDGEWID, n_var)

    # Initialise progress bar for plotting
    if print_on:
        pb.printProgressBar(0, len(x), prefix='BEGINNING PLOTTING', suffix='COMPLETE', length=40)

    # Run through x and y sets and create plots
    for i in range(len(x)):

        # Update progress bar
        if print_on:
            pb.printProgressBar(i, len(x), prefix='PLOTTING VARIABLE %d' % (i + 1), suffix='COMPLETE', length=40)

        # Plot x and y within parameters
        PLOT = plt.errorbar(np.array(x[i]), np.array(y[i]), xerr=x_error[i], yerr=y_error[i], color=COLOURS[i],
                            ms=SIZES[i], fmt=POINTSTYLES[i], mec=EDGECOLOURS[i], lw=LWID[i], elinewidth=ERBWID[i],
                            mew=EDGEWID[i])

        # Appends the artist label for this plot to HANDLES for legend entry
        if DATALABELS[i] is not None:
            HANDLES.append(PLOT[0])

        # Fills underneath line/scatters
        # WARNING: THIS IS IN ALPHA DEVELOPMENT!
        if FILL_COLOURS[i] is not None:
            plt.fill(np.array(x[i]), np.array(y[i]), alpha=SHADE[i], color=FILL_COLOURS[i])

    # Finishes progress bar
    if print_on:
        pb.printProgressBar(len(x), len(x), prefix='FINISHED PLOTTING!', suffix='COMPLETE!', length=40)

    return HANDLES


def determine_axis(x, y):
    """ Determines the min and max of x and y axis from data bounds

    Args:
        x ([[float]]): All x-axis data
        y ([[float]]): All y-axis data

    Return:
        axis_range (np.array): The dimensions of the axes

    """

    # Tells Python to edit these variables globally not locally 
    # i.e everywhere, not just in this method 
    global x_min, x_max, y_min, y_max, axis_range

    # Arrays to hold the min/max of each variable of DATASET's x and y 
    x_mins = []
    x_maxs = []
    y_mins = []
    y_maxs = []

    # Cycles through x and y to find min/max of each variable
    for i in range(len(x)):
        x_mins.append(np.min(x[i]))
        x_maxs.append(np.max(x[i]))
    for j in range(len(y)):
        y_mins.append(np.min(y[j]))
        y_maxs.append(np.max(y[j]))

    # Adds the min/max of the min/maxs to axis range
    axis_range = [np.min(x_mins), np.max(x_maxs), np.min(y_mins), np.max(y_maxs)]

    # Adds a 10% border to axis dimensions so data can be more clearly seen
    axis_range = add_borders(axis_range)

    return axis_range


def add_borders(axis_range):
    """Adds a 10% extra border to the axis dimensions

    Args:
        axis_range (np.array): Length 4 array of x and y min and maxs of axes

    Returns:
        axis_range (np.array): With a 10% extra range
    """

    for i in range(len(axis_range)):

        if i == 0 or i == 2:
            if axis_range[i] > 0:
                axis_range[i] = 0.9 * axis_range[i]
            if axis_range[i] < 0:
                axis_range[i] = 1.1 * axis_range[i]

        if i == 1 or i == 3:
            if axis_range[i] < 0:
                axis_range[i] = 0.9 * axis_range[i]
            if axis_range[i] > 0:
                axis_range[i] = 1.1 * axis_range[i]

    return axis_range


def create_figure(x, y, x_error=[None], y_error=[None], DATALABELS=[None], COLOURS=[None], FILL_COLOURS=[None],
                  SHADE=[0.5], SIZES=[None], POINTSTYLES=['none'], EDGECOLOURS=[None], EDGEWID=[2], LWID=[2],
                  ERBWID=[2], figure_name="Fig1.pdf", x_label='x-label', y_label='y-label', figsize=(10, 10),
                  axis_range=[0.0, 10.0, 0.0, 10.0], x_ticks=None, x_tick_labels=None,
                  y_ticks=None, y_tick_labels=None):
    """ Method capable of plotting multiple sets of 2D variables

    Args:
        x ([[float]]): All arrays of x-axis data to plot
        y ([[float]]): All arrays of y-axis data to plot
        x_error ([[float]]/[[tuple]]): Errors for each x value. Supports tuples for asymmetric errors
        y_error ([[float]]/[[tuple]]): Errors for each y value. Supports tuples for asymmetric errors
        DATALABELS ([str]): Labels for legend of plots
        COLOURS ([str]): Colours for each plot
        FILL_COLOURS ([str]): Colour of fill for plot. If None is entered then no fill is used
        SHADE ([float]): Alpha value for fill under each plot
        SIZES ([float]): Size of points for each plot
        POINTSTYLES ([str]): Style of point/ line for each plot
        EDGECOLOURS ([str]): Colour for points of each plot
        EDGEWID ([float]): Width of point edge for each plot
        LWID ([float]): Width of line for each plot
        ERBWID ([float]): Width of error bars for each plot

    """

    # WELCOME MESSAGE ========================================================
    if print_on:
        print("\nWELCOME TO Plot2D")
        print("PART OF THE LancAstro PACKAGE \n")

    # Sets up the figure size from parameters before plotting commences  
    plt.figure(figsize=figsize)

    # Sets the x,y-axis scales from parameters
    plt.xscale(x_scale)
    plt.yscale(y_scale)

    # Calls plot method to plot each variable from the dataset    
    HANDLES = plot(x, y, x_error=x_error, y_error=y_error, DATALABELS=DATALABELS, COLOURS=COLOURS,
                   FILL_COLOURS=FILL_COLOURS, SHADE=SHADE, SIZES=SIZES, POINTSTYLES=POINTSTYLES, EDGEWID=EDGEWID,
                   EDGECOLOURS=EDGECOLOURS, LWID=LWID, ERBWID=ERBWID)

    # If set to, automatically sets the min and max of the axis
    # from the range of all the data
    if auto_axis:
        axis_range = determine_axis(x, y)

    # Places grid if set to do so
    plt.grid(grid)

    # Adds minor ticks to axis if set to do so
    if minorticks:
        plt.minorticks_on()

    # Creates figure axis from defined parameters
    plt.axis(axis_range)

    # Creates x and y axis labels from settings 
    plt.xlabel(r'%s' % x_label, {'color': "%s" % x_colour, 'fontsize': x_size})
    plt.ylabel(r'%s' % y_label, {'color': '%s' % y_colour, 'fontsize': y_size})

    # If true, creates a legend
    if legend_on:
        if print_on:
            print("Producing legend")

        leg = plt.legend(handles=HANDLES, labels=DATALABELS)

        # Sets legend frame to be transparent if set to do so
        if not frame:
            leg.get_frame().set_facecolor('none')

        # Removes legend border if set to False
        if not frame_border:
            leg.get_frame().set_linewidth(0.0)

        # Sets frame border to defiend
        if frame_border:
            leg.get_frame().set_edgecolor(frame_colour)

    # Sets ticks on the top and right axis if true
    if both_axis:
        plt.tick_params(which="both", direction="in", top=True, right=True)
    if print_on:
        print("Saving figure to %s" % figure_name)

    if save_fig:
        plt.savefig(figure_name)

        if print_on:
            print("Figure plotted! SCIENCE!")

    if display_figure:
        plt.show()
    else:
        plt.close()


# =====================================================================================================================
#                                                     PARAMETERS
# =====================================================================================================================
# FIGURE AND AXIS SIZES ===============================================================================================
mpl.rcParams['xtick.labelsize'] = 14  # Size of x-tick labels
mpl.rcParams['ytick.labelsize'] = 14  # Size of y-tick labels
mpl.rcParams['axes.labelsize'] = 14  # Size of axes labels
mpl.rcParams['image.origin'] = 'lower'  #

mpl.rcParams['axes.linewidth'] = 2.5  # Size of axes line width
mpl.rcParams['xtick.major.size'] = 14  # Size of major x-ticks
mpl.rcParams['xtick.minor.size'] = 6  # Size of minor x-ticks
mpl.rcParams['xtick.major.width'] = 2.5  # Width of major x-ticks
mpl.rcParams['xtick.minor.width'] = 1.5  # Width of minor x-ticks
mpl.rcParams['ytick.major.size'] = 14  # Size of major y-ticks
mpl.rcParams['ytick.minor.size'] = 6  # Size of minor y-ticks
mpl.rcParams['ytick.major.width'] = 2.5  # Width of major y-ticks
mpl.rcParams['ytick.minor.width'] = 1.5  # Width of minor y-ticks

# Sets font style and size
mpl.rcParams.update({'font.size': 18, 'font.weight': 'bold'})

# FIGURE PARAMETERS ===================================================================================================
share_plot = True  # If true, plot all variables onto same axis. If false, outputs each plot in separate figures
figure_name = 'Test.pdf'  # Defines default filename

figsize = (10, 7)  # Defines figure size

display_figure = True  # Opens figure if true

save_fig = True  # Saves the figure to file if true

print_on = False  # Includes print statements in output if true

# POINTS & FILL PARAMETERS ============================================================================================
POINTSTYLES = ['*', 'o']  # Point style
SIZES = [1, 2]  # Point sizes
COLOURS = ['0.3', 'b']  # Point fill colours
EDGECOLOURS = ['face', 'g']  # Point edge colours

# AXIS PARAMETERS =====================================================================================================
x_min = 0.0  # Minimum x-axis value to plot
x_max = 10.0  # Maximum x-axis value to plot
y_min = 0.0  # Minimum y-axis value to plot
y_max = 10.0  # Maximum y-axis value to plot

axis_range = [x_min, x_max, y_min, y_max]

x_scale = 'linear'  # Defines x-axis scaling
y_scale = 'linear'  # Defines y-axis scaling

auto_axis = True  # Auto define dimesions of axises

both_axis = True  # Places ticks on top and right axes too

auto_tick = True  # Auto create the ticks
# AXIS LABELS =========================================================================================================
x_label = "x-label"  # x-axis label
y_label = "y-label"  # y-axis label
x_size = 16  # Size of x-axis label
y_size = 16  # Size of y-axis label
x_colour = "k"  # Colour of x-axis label
y_colour = "k"  # Colour of y-axis label

minorticks = True  # Places minor ticks on axis if True
grid = False  # Places grid on plot if True

# LEGEND PARAMETERS ===================================================================================================
legend_on = True  # Creates legend if true

frame_border = False  # Puts border around the legend
frame = True  # Sets frame to visible
frame_colour = "k"  # Sets frame colour

# END PROGRAM =========================================================================================================
