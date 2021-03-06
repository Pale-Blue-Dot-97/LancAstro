# =====================================================================================================================
# =================================================== MULTIFIG ========================================================
# =====================================================================================================================
"""Script for creating a grid of subplots using Plot2D.py

Author: Harry Baker
Date: 02.10.2019
Version: 1.0
== RELEASE ==

Part of the LancAstro.py project for XGAL

Example:
    $ MFTest.py

"""

# =====================================================================================================================
#                                                     IMPORTS
# =====================================================================================================================

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import Dependencies.progress_bar as pb
import Plot2D as laplt


# =====================================================================================================================
#                                                     FUNCTIONS
# =====================================================================================================================


def create_grid(x, y, shape, x_error=[[None]], y_error=[[None]], DATALABELS=[[None]], COLOURS=[[None]],
                FILL_COLOURS=[[None]], SHADE=[[0.5]], SIZES=[[None]], POINTSTYLES=[['none']], EDGECOLOURS=[[None]],
                EDGEWID=[[2]], LWID=[[2]], ERBWID=[[2]], figure_name="MultiFig.pdf", x_label='x-label',
                y_label='y-label', figsize=(10, 10), axis_range=[0.0, 10.0, 0.0, 10.0], x_ticks=None,
                x_tick_label=None, y_ticks=None, y_tick_label=None, num_x_tick=10, num_y_tick=10):
    """ Creates a figure (or multiple figures) of the plots of data supplied

    Args:
        x ([[[float]]]): Data to plot on x-axis
        y ([[[float]]]): Data to plot on y-axis
        shape ([[float]]): Grid defining the positioning of each subplot
        x_error ([[[float]]]/[[[tuple]]]): Errors for each x value. Supports tuples for asymmetric errors
        y_error ([[[float]]]/[[[tuple]]]): Errors for each y value. Supports tuples for asymmetric errors
        DATALABELS ([[str]]): Label for each variable for legend
        COLOURS ([[str]]): Colours for each plot
        FILL_COLOURS ([[str]]): Colour of fill for plot. If None is entered then no fill is used
        SHADE ([[float]]): Alpha value for fill under each plot
        SIZES ([[float]]): Size of points for each plot
        POINTSTYLES ([[str]]): Style of point/ line for each plot
        EDGECOLOURS ([[str]]): Colour for points of each plot
        EDGEWID ([[float]]): Width of point edge for each plot
        LWID ([[float]]): Width of line for each plot
        ERBWID ([[float]]): Width of error bars for each plot
        figure_name (str): Name of figure
        x_label (str): x-axis label
        y_label (str): y_axis label
        figsize (tuple): Figure size
        axis_range (np.array): Dimensions of subplot axises
        x_ticks (list/array/np.array): Positions of x-axis ticks
        x_tick_label ([str]): Labels for each x-tick
        y_ticks (list/array/np.array): Positions of y-axis ticks
        y_tick_label ([str]): Labels for each y-tick
        num_x_tick (int): Number of desired x-axis ticks
        num_y_tick (int): Number of desired y-axis ticks

    """

    # WELCOME MESSAGE
    if print_on:
        print("\nWELCOME TO MULTIFIG")
        print("PART OF THE LancAstro PACKAGE \n")

    # Sets up the figure size from parameters before plotting commences  
    plt.figure(figsize=figsize)

    # Sets the x,y-axis scales from parameters
    plt.xscale(x_scale)
    plt.yscale(y_scale)

    n_plots = len(x)

    n_cols = len(shape)
    n_rows = len(shape[0])

    x_error = laplt.default_lengths(x_error, n_plots)
    y_error = laplt.default_lengths(y_error, n_plots)
    DATALABELS = laplt.default_lengths(DATALABELS, n_plots)
    COLOURS = laplt.default_lengths(COLOURS, n_plots)
    FILL_COLOURS = laplt.default_lengths(FILL_COLOURS, n_plots)
    SIZES = laplt.default_lengths(SIZES, n_plots)
    POINTSTYLES = laplt.default_lengths(POINTSTYLES, n_plots)
    EDGECOLOURS = laplt.default_lengths(EDGECOLOURS, n_plots)
    SHADE = laplt.default_lengths(SHADE, n_plots)
    LWID = laplt.default_lengths(LWID, n_plots)
    ERBWID = laplt.default_lengths(ERBWID, n_plots)
    EDGEWID = laplt.default_lengths(EDGEWID, n_plots)

    left_edge = []
    bottom_edge = []

    for j in range(len(shape)):
        left_edge.append(shape[j][0])

    for j in np.arange(0, n_rows, 1):
        bottom_edge.append(shape[len(shape) - 1][j])

    for i in range(len(x)):

        plt.subplot(n_cols, n_rows, i + 1)

        # Calls plot method to plot each variable from the dataset    
        HANDLES = laplt.plot(x[i], y[i], x_error=x_error[i], y_error=y_error[i], DATALABELS=DATALABELS[i],
                             COLOURS=COLOURS[i], FILL_COLOURS=FILL_COLOURS[i], SHADE=SHADE[i], SIZES=SIZES[i],
                             POINTSTYLES=POINTSTYLES[i], EDGEWID=EDGEWID[i], EDGECOLOURS=EDGECOLOURS[i], LWID=LWID[i],
                             ERBWID=ERBWID[i])

        # If set to, automatically sets the min and max of the axis
        # from the range of all the data
        if auto_axis:
            axis_range = laplt.determine_axis(x[i], y[i])

        if auto_tick:
            x_ticks, y_ticks = laplt.construct_ticks(axis_range, num_x_tick, num_y_tick)

        # Places grid if set to do so
        plt.grid(grid)

        # Adds minor ticks to axis if set to do so
        if minorticks:
            plt.minorticks_on()

        # Creates figure axis from defined parameters
        plt.axis(axis_range)

        if i + 1 in left_edge:
            plt.ylabel(r'%s' % y_label, {'color': '%s' % y_colour, 'fontsize': y_size})
            plt.yticks(y_ticks, y_tick_label)

        if i + 1 in bottom_edge:
            plt.xlabel(r'%s' % x_label, {'color': "%s" % x_colour, 'fontsize': x_size})
            plt.xticks(x_ticks, x_tick_label)

        if i + 1 not in left_edge:
            plt.yticks(y_ticks, [])

        if i + 1 not in bottom_edge:
            plt.xticks(x_ticks, [])

        # If true, creates a legend
        if legend_on:
            if print_on:
                print("Producing legend")

            leg = plt.legend(handles=HANDLES, labels=DATALABELS[i])

            # Sets legend frame to be transparent if set to do so
            if not frame:
                leg.get_frame().set_facecolor('none')

            # Removes legend border if set to False
            if not frame_border:
                leg.get_frame().set_linewidth(0.0)

            # Sets frame border to defined colour
            if frame_border:
                leg.get_frame().set_edgecolor(frame_colour)

        # Sets ticks on the top and right axis if true
        if both_axis:
            plt.tick_params(which="both", direction="in", top=True, right=True)

    # Adjust plots so they don't have any white space in between them.
    plt.subplots_adjust(wspace=0, hspace=0)

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
#                                                   PARAMETERS
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
mpl.rcParams.update({'font.size': 10, 'font.weight': 'bold'})

# FIGURE PARAMETERS ===================================================================================================
share_plot = True  # If true, plot all variables onto same axis
# if false, outputs each plot in seperate
# figures
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
# x_min = 0.0                 # Minimum x-axis value to plot
# x_max = 10.0                # Maximum x-axis value to plot
# y_min = 0.0                 # Minimum y-axis value to plot
# y_max = 10.0                # Maximum y-axis value to plot

# axis_range = [x_min,x_max,y_min,y_max]

x_scale = 'linear'  # Defines x-axis scaling
y_scale = 'linear'  # Defines y-axis scaling

auto_axis = True  # Auto define dimesions of axises

auto_tick = True # Automatically calculates positions of axis ticks

both_axis = True  # Places ticks on top and right axes too

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
