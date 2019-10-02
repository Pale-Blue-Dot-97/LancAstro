# =====================================================================================================================
# ===================================================== PLOT3D ========================================================
# =====================================================================================================================
"""
Script for 3D plotting of data stored in FITS files, ascii tables

Author: Harry Baker
Date: 02.10.2019
Version: 0.8.1
== BETA ==

Part of the LancAstro.py project for XGAL

Example:


"""

# =====================================================================================================================
#                                                     IMPORTS
# =====================================================================================================================

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import random
import Dependencies.progress_bar as pb
import Plot2D as laplt2D


# =====================================================================================================================
#                                                    FUNCTIONS
# =====================================================================================================================

def construct_ticks(axis_range, num_x, num_y):
    x_step = (axis_range[1] - axis_range[0]) / float(num_x)
    y_step = (axis_range[3] - axis_range[2]) / float(num_y)

    x_ticks = np.round(np.arange(axis_range[0], axis_range[1], x_step), 1)
    y_ticks = np.round(np.arange(axis_range[2], axis_range[3], y_step), 1)

    return x_ticks, y_ticks


def plot(ax, x, y, z, DATALABELS=[None], COLOURS=[None], SIZES=[2], POINTSTYLES=['none'], CMAP=[None]):
    """ Method capable of plotting multiple sets of 2D variables

    Args:
        x ([[float]]): All arrays of x-axis data to plot
        y ([[float]]): All arrays of y-axis data to plot
        DATALABELS ([str]): Labels for legend of plots
        COLOURS ([str]): Colours for each plot
        SIZES ([float]): Size of points for each plot
        POINTSTYLES ([str]): Style of point/ line for each plot
        CMAP ():

    Returns:
        HANDLES: Handles for each plot to be used to construct a legend and final figure

    """

    HANDLES = []

    if not isinstance(x[0], (list, np.ndarray)): x = [x]
    if not isinstance(y[0], (list, np.ndarray)): y = [y]
    if not isinstance(z[0], (list, np.ndarray)): z = [z]

    n_var = len(x)

    # Makes sure all parameters are the same length as number of plots
    DATALABELS = laplt2D.default_lengths(DATALABELS, n_var)
    COLOURS = laplt2D.default_lengths(COLOURS, n_var)
    SIZES = laplt2D.default_lengths(SIZES, n_var)
    POINTSTYLES = laplt2D.default_lengths(POINTSTYLES, n_var)

    # Initialise progress bar for plotting
    if print_on:
        pb.printProgressBar(0, len(x), prefix='BEGINNING PLOTTING', suffix='COMPLETE', length=40)

    # Run through x and y sets and create plots
    for i in range(len(x)):

        # Update progress bar
        if print_on:
            pb.printProgressBar(i, len(x), prefix='PLOTTING VARIABLE %d' % (i + 1), suffix='COMPLETE', length=40)
        c = None

        if COLOURS[i] == 'colourbar' and CMAP[i] is not None:
            c = np.array(z[i])

        # Plot x and y within parameters
        PLOT = ax.scatter(xs=np.array(x[i]), ys=np.array(y[i]), zs=np.array(z[i]), c=c, cmap=CMAP[i], s=SIZES[i],
                          marker=POINTSTYLES[i])

        cb = plt.colorbar(PLOT)

        if not visible_colourbar:
            cb.remove()

        # Appends the artist label for this plot to HANDLES for legend entry
        if DATALABELS[i] is not None:
            HANDLES.append(PLOT)

    # Finishes progress bar
    if print_on:
        pb.printProgressBar(len(x), len(x), prefix='FINISHED PLOTTING!', suffix='COMPLETE!', length=40)
        print("\n")

    return ax, HANDLES


def determine_axis(x, y, z=None, w=0.1):
    """ Determines the min and max of x and y axis from data bounds

    Args:
        x: All x-axis data
        y: All y-axis data
        z: All z-axis data

    Return:
        axis_range: The dimensions of the axes
    """

    # Tells Python to edit these variables globally not locally 
    # i.e everywhere, not just in this method 
    # global x_min,x_max,y_min,y_max,axis_range

    # Arrays to hold the min/max of each variable of DATASET's x and y 
    x_mins = []
    x_maxs = []
    y_mins = []
    y_maxs = []
    z_mins = []
    z_maxs = []

    # Cycles through x and y to find min/max of each variable
    for i in range(len(x)):
        x_mins.append(np.min(x[i]))
        x_maxs.append(np.max(x[i]))
    for j in range(len(y)):
        y_mins.append(np.min(y[j]))
        y_maxs.append(np.max(y[j]))

    if z is not None:
        for k in range(len(z)):
            z_mins.append(np.min(z[k]))
            z_maxs.append(np.max(z[k]))

    # Adds the min/max of the min/maxs to axis range
    axis_range = [np.min(x_mins), np.max(x_maxs), np.min(y_mins), np.max(y_maxs),
                  np.min(z_mins), np.max(z_maxs)]

    # Adds a 10% border to axis dimensions so data can be more clearly seen
    axis_range = add_borders(axis_range, w)

    return axis_range


def add_borders(axis_range, w):
    """ Adds a 10% extra border to the axis dimensions

    Args:
        axis_range: Length 4 array of x and y min and maxs of axes

    Returns:
        axis_range: With a 10% extra range 
    """

    for i in range(len(axis_range)):
        try:
            if i == 0 or i == 2 or i == 4:
                if axis_range[i] > 0.0:
                    axis_range[i] = (1.0 - w) * axis_range[i]
                if axis_range[i] < 0.0:
                    axis_range[i] = (1.0 + w) * axis_range[i]

            if i == 1 or i == 3 or i == 5:
                if axis_range[i] < 0.0:
                    axis_range[i] = (1.0 - w) * axis_range[i]
                if axis_range[i] > 0.0:
                    axis_range[i] = (1.0 + w) * axis_range[i]
        except:
            continue

    return axis_range


def create_figure(x, y, z, DATALABELS=[None], COLOURS=[None], CMAP=[None], SIZES=[2], POINTSTYLES=['none'],
                  figure_name="Plot3D_fig.pdf", x_label=None, y_label=None, z_label=None, figsize=(10., 10.),
                  axis_range=[0.0, 10.0, 0.0, 10.0], label_rot=[20.0, 340.0, 90.0], elev=10., azim=60.):
    """ Creates a figure (or multiple figures) of the plots of data supplied

    Args:
        x ([[float]]): All arrays of x-axis data to plot
        y ([[float]]): All arrays of y-axis data to plot
        z ([[float]]): All arrays of z-axis data to plot
        DATALABELS ([str]): Labels for legend of plots
        COLOURS ([str]): Colours for each plot

        SIZES ([float]): Size of points for each plot
        POINTSTYLES ([str]): Style of point/ line for each plot

    """

    # WELCOME MESSAGE ========================================================
    if print_on:
        print("\nWELCOME TO Plot3D")
        print("PART OF THE LancAstro PACKAGE \n")

    # Sets up the figure size from parameters before plotting commences  
    fig = plt.figure(figsize=figsize)

    ax = fig.add_subplot(111, projection='3d')

    # Sets the x,y-axis scales from parameters
    ax.set_xscale(x_scale)
    ax.set_yscale(y_scale)
    ax.set_zscale(z_scale)

    # Calls plot method to plot each variable from the dataset    
    ax, HANDLES = plot(ax, x, y, z, DATALABELS=DATALABELS, COLOURS=COLOURS, SIZES=SIZES, POINTSTYLES=POINTSTYLES,
                       CMAP=CMAP)

    # If set to, automatically sets the min and max of the axis
    # from the range of all the data
    if auto_axis:
        axis_range = determine_axis(x, y, z=z, w=w)

    # Sets elevation and azimuthal angle of plot 
    ax.view_init(elev=elev, azim=azim)

    # Places grid if set to do so
    ax.grid(grid)

    # Creates figure axis from defined parameters
    ax.set_xlim(axis_range[0], axis_range[1])
    ax.set_ylim(axis_range[2], axis_range[3])
    ax.set_zlim(axis_range[4], axis_range[5])

    # Creates x and y axis labels from settings 
    if visible_axes:
        ax.set_xlabel(x_label, color=x_colour, fontsize=x_size, rotation=label_rot[0])
        ax.set_ylabel(y_label, color=y_colour, fontsize=y_size, rotation=label_rot[1])
        ax.set_zlabel(z_label, color=z_colour, fontsize=z_size, rotation=label_rot[2])

    if not visible_axes:
        fig.axes.get_xaxis().set_visible(False)
        fig.axes.get_yaxis().set_visible(False)
        fig.axes.get_zaxis().set_visible(False)

    # If true, creates a legend
    if legend_on:
        if print_on:
            print("Producing legend")

        leg = fig.legend(handles=HANDLES, labels=DATALABELS)

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
        ax.tick_params(which="both", direction="in", top=True, right=True)

    if print_on:
        print("Saving figure to %s" % figure_name)

    if save_fig:
        fig.savefig(figure_name, bbox_inches='tight')

        if print_on:
            print("Figure plotted! SCIENCE!")

    if display_figure:
        plt.show()
    if display_figure:
        plt.close()


def create_visualisation(x, y, z, DATALABELS=[None], COLOURS=[None], CMAP=[None], SIZES=[2], POINTSTYLES=['none'],
                         x_label=None, y_label=None, z_label=None, figure_name="Plot3D_fig.pdf", figsize=(10., 10.),
                         axis_range=[0.0, 10.0, 0.0, 10.0], label_rot=[20.0, 340.0, 90.0], elev=None, azim=None,
                         no_frames=20, folder='Frames'):
    """ Creates a figure (or multiple figures) of the plots of data supplied

    Args:
        x (2D Array of floats): Data to plot on x-axis

        y (2D Array of floats): Data to plot on y-axis

        x_error (2D Array of floats/tuples): Errors for each x value 
        Supports tuples for assymetric errors 

        y_error (2D Array of floats/tuples): Errors for each y value 
        Supports tuples for assymetric errors

        DATALABELS (Array of strings): Label for each variable for legend

        COLOURS (Array of strings): Colours    
    """

    # WELCOME MESSAGE ========================================================
    if print_on == True:
        print("\nWELCOME TO Plot3D")
        print("PART OF THE LancAstro PACKAGE \n")

    # Sets up the figure size from parameters before plotting commences  
    fig = plt.figure(figsize=figsize)

    ax = fig.add_subplot(111, projection='3d')

    # Sets the x,y-axis scales from parameters
    ax.set_xscale(x_scale)
    ax.set_yscale(y_scale)
    ax.set_zscale(z_scale)

    # Calls plot method to plot each variable from the dataset    
    ax, HANDLES = plot(ax, x, y, z, DATALABELS=DATALABELS, COLOURS=COLOURS,
                       SIZES=SIZES, POINTSTYLES=POINTSTYLES, CMAP=CMAP)

    # If set to, automatically sets the min and max of the axis
    # from the range of all the data
    if auto_axis == True:
        axis_range = determine_axis(x, y, z=z, w=w)

    # Places grid if set to do so
    ax.grid(grid)

    ax.axis(visible_axes)

    # Creates figure axis from defined parameters
    ax.set_xlim(axis_range[0], axis_range[1])
    ax.set_ylim(axis_range[2], axis_range[3])
    ax.set_zlim(axis_range[4], axis_range[5])

    # Creates x and y axis labels from settings 
    if visible_axes == True:
        ax.set_xlabel(x_label, color=x_colour, fontsize=x_size, rotation=label_rot[0])
        ax.set_ylabel(y_label, color=y_colour, fontsize=y_size, rotation=label_rot[1])
        ax.set_zlabel(z_label, color=z_colour, fontsize=z_size, rotation=label_rot[2])

    if visible_axes == False:
        ax.set_xlabel([])
        ax.set_ylabel([])
        ax.set_zlabel([])
        # ax.axes.get_xaxis().set_visible(False)
        # ax.axes.get_yaxis().set_visible(False)
        # ax.axes.get_zaxis().set_visible(False)

    # If true, creates a legend
    if legend_on == True:
        if print_on == True:
            print("Producing legend")

        leg = fig.legend(handles=HANDLES, labels=DATALABELS)

        # Sets legend frame to be transparent if set to do so
        if frame == False:
            leg.get_frame().set_facecolor('none')

        # Removes legend border if set to False
        if frame_border == False:
            leg.get_frame().set_linewidth(0.0)

        # Sets frame border to defiend
        if frame_border == True:
            leg.get_frame().set_edgecolor(frame_colour)

    # Sets ticks on the top and right axis if true
    if both_axis == True:
        ax.tick_params(which="both", direction="in", top=True, right=True)

    # Initialise progress bar for plotting
    if print_on == True:
        pb.printProgressBar(0, no_frames + 2, prefix='BEGINNING VISUALISATION',
                            suffix='COMPLETE', length=40)

    for i in np.arange(0, no_frames, 1):

        cpfig = fig

        # Sets elevation and azimuthal angle of plot
        if azim == None:
            ax.view_init(elev=elev, azim=(float(i) / float(no_frames) * 360.0))
        if elev == None:
            ax.view_init(elev=(float(i) / float(no_frames) * 360.0), azim=azim)
        if elev == None and azim == None:
            ax.view_init(elev=(float(i) / float(no_frames) * 360.0), azim=(float(i) / float(no_frames) * 360.0))

        cpfig.savefig("%s%d%s" % (folder, i, figure_name), bbox_inches='tight')

        if print_on == True:
            pb.printProgressBar(i + 1, no_frames + 2, prefix='FRAME NO.%d SAVED      ' % i,
                                suffix='COMPLETE', length=40)

        plt.close()

    # Finishes progress bar
    if print_on == True:
        pb.printProgressBar(no_frames + 2, no_frames + 2, prefix='VISUALISATION COMPLETE!',
                            suffix='COMPLETE!', length=40)
        print("\n")

    # ============================================================================


#                           PARAMETERS
# ============================================================================
# FIGURE AND AXIS SIZES ======================================================
mpl.rcParams['xtick.labelsize'] = 10  # Size of x-tick labels
mpl.rcParams['ytick.labelsize'] = 10  # Size of y-tick labels

mpl.rcParams['axes.linewidth'] = 1  # Size of axes line width
mpl.rcParams['axes.labelsize'] = 10  # Size of axes labels
mpl.rcParams['image.origin'] = 'lower'  #

mpl.rcParams['xtick.major.size'] = 10  # Size of major x-ticks
mpl.rcParams['xtick.minor.size'] = 5  # Size of minor x-ticks
mpl.rcParams['xtick.major.width'] = 1.5  # Width of major x-ticks
mpl.rcParams['xtick.minor.width'] = 1  # Width of minor x-ticks

mpl.rcParams['ytick.major.size'] = 10  # Size of major y-ticks
mpl.rcParams['ytick.minor.size'] = 5  # Size of minor y-ticks
mpl.rcParams['ytick.major.width'] = 1.5  # Width of major y-ticks
mpl.rcParams['ytick.minor.width'] = 1  # Width of minor y-ticks

# Sets font style and size
mpl.rcParams.update({'font.size': 14, 'font.weight': 'bold'})

# FIGURE PARAMETERS ==========================================================
share_plot = True  # If true, plot all variables onto same axis
# if false, outputs each plot in seperate
# figures

display_figure = True  # Opens figure if true

save_fig = True  # Saves the figure to file if true

print_on = True  # Includes print statements in output if true

# AXIS PARAMETERS ============================================================
x_scale = 'linear'  # Defines x-axis scaling
y_scale = 'linear'  # Defines y-axis scaling
z_scale = 'linear'  # Defines y-axis scaling

auto_axis = True  # Auto define dimesions of axises

both_axis = True  # Places ticks on top and right axes too

auto_tick = True  # Auto create the ticks

w = 0.0  # Width of border of axis-range
# as a fraction of axis-range

visible_axes = False  # Sets visibility of axes

visible_colourbar = True

# AXIS LABELS ================================================================
x_size = 16  # Size of x-axis label
y_size = 16  # Size of y-axis label
z_size = 16  # Size of z-axis label
x_colour = "k"  # Colour of x-axis label
y_colour = "k"  # Colour of y-axis label
z_colour = "k"  # Colour of z-axis label

minorticks = True  # Places minor ticks on axis if True
grid = False  # Places grid on plot if True

# LEGEND PARAMETERS ==========================================================
legend_on = False  # Creates legend if true

frame_border = False  # Puts border around the legend
frame = True  # Sets frame to visible
frame_colour = "k"  # Sets frame colour

# END PROGRAM ================================================================
