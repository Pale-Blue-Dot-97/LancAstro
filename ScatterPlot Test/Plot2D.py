#=============================================================================
#================================ PLOT2D =====================================
#=============================================================================
#
# Script for 2D plotting of data stored in FITS files, ascii tables
#
# Author: Harry Baker
# Date: 24.07.2019
# Version: 1.1
# == RELEASE ==
#
# Part of the LancAstro.py project for XGAL
# 
# CHANGE LOG:
#   V1.1:
#     -Added convert_to_array method to check and correct if supplied
#      parameters are arrays
#     -Added save_fig option to parameters to toggle saving of figures to file
#     -Added SHADE parameter array so plots can have fills underneath
#   
#   V1.0:
#     -Renamed ScatterPlot to Plot2D for its Version 1 release for LancAstro!
#     -Added print_on option to parameters to toggle print to terminal
#     -Added display_figure option to parameters to toggle the figure pop-up
#     -Added x_label and y_label option args to create_figure method so that
#      user can define axis labels from invoking of method
#     -Added axis_range optional args to create_figure method so that user
#      can define axis dimensions from invoking of method
#     -Added legend_on option to parameters to toggle creation of legend
#     -Fixed issue with plot method where a None arguememt for 'fmt' option
#      of pyplot.errorbar is depreciated
#     -Converted plot method to support scatter, errorbars and line plots   
#
version = "1.1"
# ============================================================================
#                             IMPORTS   
# ============================================================================

import os, sys
import matplotlib as mpl
import matplotlib.pyplot as plt
from astropy.io import fits as pyfits
import numpy as np 
import random
import Dependencies.ascii_read as ar
import Dependencies.progress_bar as pb


############################ DATA LOADING ####################################
def data_load(DATANAME, COLUMNNAME,PATH=['']):
    """Load data in from columns in a FITS file
    
    Args:
        DATANAME: Array of strings with names of the FITS file containing the 
                data to be loaded
        COLUMNNAME: 2D array of strings with first axis being of equal length
                  to DATANAME, and the second axis being the names of the
                  columns to be loaded from each FITS file in DATANAME
        PATH: Optional variable that is a list of strings for paths to each
              FITS file. Should be same length as DATANAME

    Returns:
        A 2D array of data requested from FITS files

    """
    # Defines the array that will hold all the variables to plot
    DATASET = []

    # Checks that PATH is the same length as DATANAME and corrects if not
    PATH = default_lengths(PATH,len(DATANAME))

    if print_on == True:
        # Intialises the progress bar for loading the data
        pb.printProgressBar(0,2*len(DATANAME)*sum(len(x) for x in COLUMNNAME), 
                           prefix='LOADING DATA:',suffix='COMPLETE',length=40)
        
    # Loops through the FITS files
    for i in range(len(DATANAME)):
        # Exception handling in case file cannot be opened
        try:
            # Opens file into memory
            FILE = pyfits.open("%s%s"%(PATH[i],DATANAME[i]))

            if print_on == True:
                # Updates progress bar
                pb.printProgressBar(
                    i,2*len(DATANAME)*sum(len(x) for x in COLUMNNAME), 
                    prefix = 'OPENING DATA FILE %s%s:'%(PATH[i],DATANAME[i]), 
                    suffix = 'COMPLETE',length=40)

            # Loads Extension 1 of the FITS file into memory
            TABLE = FILE[1].data

            # Loops through and loads all the columns in the file to be loaded
            for j in range(len(COLUMNNAME[i])):

                if print_on == True:
                    # Progress bar updated to reflect this
                    pb.printProgressBar(
                        i+j,2*len(DATANAME)*sum(len(x) for x in COLUMNNAME), 
                        prefix = 'LOADING COLUMN %s:'%COLUMNNAME[i][j], 
                        suffix = 'COMPLETE',length=40)

                # Loads column into memory
                VAR = TABLE.field("%s"%COLUMNNAME[i][j])

                # Removes any 'nan' (Not an Actual Number) strings that are
                # placeholders if there are blank entries as these would throw
                # errors
                CleanVAR = [x for x in VAR if str(x) != 'nan']

                # Adds 'cleaned' variable into DATASET
                DATASET.append(CleanVAR) 

        # If file cannot be opened, throws an exception   
        except:
            print("FILE COULD NOT BE OPENED")

    if print_on == True:
        # Completes progress bar once all variables are loaded into DATASET
        pb.printProgressBar(2*len(DATANAME)*sum(len(x) for x in COLUMNNAME),
                            2*len(DATANAME)*sum(len(x) for x in COLUMNNAME), 
                            prefix = 'DATASET LOADED READY FOR PLOTTING!', 
                            suffix = 'COMPLETE!',length=40)

    return DATASET # Returns 2D array of variables

###################### FIND MIN/MAX OF DATA ################################
def _get_outer_edges(a):
    """
    Determine the outer bin edges to use, from either the data or the range
    argument
    ########
    Edited from Numpy source code:
    https://github.com/numpy/numpy/blob/master/numpy/lib/histograms.py
    """
    
    RANGE = np.min(a), np.max(a)
    if RANGE is not None:
        first_edge, last_edge = RANGE
        if first_edge > last_edge:
            raise ValueError(
                'Maximum value must be larger than minimum in range parameter.')
        if not (np.isfinite(first_edge) and np.isfinite(last_edge)):
            raise ValueError(
                "Supplied range of [{}, {}] is not finite".format(first_edge, last_edge))
    elif a.size == 0:
        # handle empty arrays. Can't determine range, so use 0-1.
        first_edge, last_edge = 0, 1
    else:
        first_edge, last_edge = np.min(a), np.max(a)
        if not (np.isfinite(first_edge) and np.isfinite(last_edge)):
            raise ValueError(
                "Autodetected range of [{}, {}] is not finite".format(first_edge, last_edge))

    # expand empty range to avoid divide by zero
    if first_edge == last_edge:
        first_edge = first_edge - 0.5
        last_edge = last_edge + 0.5

    return first_edge, last_edge

def convert_to_array(var):
    """ Checks if parameter supplied is an array and corrects if not

    Args:
        var: parameter value/array to be checked and corrected

    Return:
        var
    """

    if isinstance(var,(list,tuple,np.ndarray)) == False:
        var = np.array(var)

    return var

def default_lengths(param, n_var):
    """
    Checks every parameter list is the same length
    i.e if user only defines min and max for 1 variable
    Corrects these to all have same length 

    Args:
        param: Parameter array to be set to correct length

        n_var: Number of variables therefore the length to set param to

    Return:
        param 
    """

    param = convert_to_array(param)

    if len(param) != n_var:
        if print_on == True:
            print("WARNING: Not every variable has defined range limits!")
            print("Setting all variables to default of 1st entry \n")
        
        for i in range(len(param)):
            param[i] = param[0]

        for i in range(n_var - len(param)):
            param.append(param[0])

    return param

def colour_cycle(n_var,COLOURS):
    if len(COLOURS) != n_var:
        for i in range(n_var):
            r = lambda: random.randint(0,255)
            random_colour = '#%02X%02X%02X' % (r(),r(),r())
            if i < len(COLOURS):
                COLOURS[i] = random_colour
            else:
                COLOURS.append(random_colour)

def plot(x,y,x_error=[None],y_error=[None],DATALABELS=[None],COLOURS=[None],
         FILL_COLOURS=[None],SHADE=[0.5],SIZES=[None],POINTSTYLES=['none'],
         EDGECOLOURS=[None],LWID=[0]):
    
    HANDLES = []

    x = convert_to_array(x)
    y = convert_to_array(y)

    n_var = len(x)

    x_error = default_lengths(x_error,n_var)
    y_error = default_lengths(y_error,n_var)
    DATALABELS = default_lengths(DATALABELS,n_var)
    COLOURS = default_lengths(COLOURS,n_var)
    FILL_COLOURS = default_lengths(FILL_COLOURS,n_var)
    SIZES = default_lengths(SIZES,n_var)
    POINTSTYLES = default_lengths(POINTSTYLES,n_var)
    EDGECOLOURS = default_lengths(EDGECOLOURS,n_var)
    SHADE = default_lengths(SHADE,n_var)
    LWID = default_lengths(LWID,n_var)

    if print_on == True:
        pb.printProgressBar(0,len(x),prefix='BEGINNING PLOTTING',
                            suffix='COMPLETE', length=40)

    for i in range(len(x)):
        if print_on == True:
            pb.printProgressBar(i,len(x),prefix='PLOTTING VARIABLE %d'%(i+1),
                                suffix='COMPLETE',length=40)
        
        PLOT = plt.errorbar(np.array(x[i]),np.array(y[i]),xerr=x_error[i],
                            yerr=y_error[i],color=COLOURS[i],ms=SIZES[i],
                            fmt=POINTSTYLES[i],edgecolor=EDGECOLOURS[i],
                            lw=LWID[i],label=None)  
        
        if DATALABELS[i] != None:
            HANDLES.append(PLOT[0])

        if FILL_COLOURS[i] != None:
            plt.fill(np.array(x[i]),np.array(y[i]),alpha=SHADE[i],
                     color=FILL_COLOURS[i])

    if print_on == True:
        pb.printProgressBar(len(x),len(x),prefix='FINISHED PLOTTING!',
                            suffix='COMPLETE!',length=40)    
        
    return HANDLES


def determine_axis(x,y):

    global x_min,x_max,y_min,y_max,axis_range

    x_mins = []
    x_maxs = []
    y_mins = []
    y_maxs = []
    
    for i in range(len(x)):
        x_mins.append(np.min(x[i]))
        x_maxs.append(np.max(x[i]))
    for j in range(len(y)):
        y_mins.append(np.min(y[j]))
        y_maxs.append(np.max(y[j]))

    if y_max < 0:
        y_max = 0.9*y_max
    if y_max > 0:
        y_max = 1.1*y_max

    axis_range = [np.min(x_mins),np.max(x_maxs),np.min(y_mins),np.max(y_maxs)]

    return axis_range

#############################################################################

def create_figure(x,y,x_error=[None],y_error=[None],DATALABELS=[None],
                  COLOURS=[None],FILL_COLOURS=[None],SHADE=[0.5],SIZES=[None],
                  POINTSTYLES=['none'],EDGECOLOURS=[None],LWID=[0],
                  figure_name="Fig1.pdf", x_label='x-label',
                  y_label='y-label',axis_range=[0.0,10.0,0.0,10.0],
                  txt_labels=[[None]]):
    """ Creates a figure (or multiple figures) of the plots of data supplied

    Args:
        x (2D Array of floats): Data to plot on x-axis

        y (2D Array of floats): Data to plot on y-axis

        x_error (2D Array of floats/tuples): Errors for each x value 
        Supports tuples for assymetric errors 

        y_error (2D Array of floats/tuples): Errors for each y value 
        Supports tuples for assymetric errors

        DATALABELS (Array of strings): Label for each variable for legend

        COLOURS (Array of strings):   
    """

    # WELCOME MESSAGE ========================================================
    if print_on == True:    
        print("\nWELCOME TO Plot2D %s"%version)
        print("PART OF THE LancAstro PACKAGE \n")


    #  
    plt.figure(figsize=figsize)
    
    plt.xscale(x_scale)
    plt.yscale(y_scale)

    # Calls plot method to plot each variable from the dataset
    # plot() calls methods to construct bins    
    HANDLES = plot(x,y,x_error,y_error,DATALABELS,COLOURS,FILL_COLOURS,SHADE,
                   SIZES,POINTSTYLES,EDGECOLOURS,LWID)

    # If set to, automatically sets the min and max of the axis
    # from the range of all the data
    if auto_axis == True:
        axis_range = determine_axis(x,y)

    # Places grid if set to do so
    plt.grid(grid)

    # Adds minor ticks to axis if set to do so
    if minorticks == True:
        plt.minorticks_on()

    # Creates figure axis from defined parameters
    plt.axis(axis_range)

    # Creates x and y axis labels from settings 
    plt.xlabel(r'%s'%x_label, {'color'    : "%s"%x_colour,  'fontsize'   : x_size })
    plt.ylabel(r'%s'%y_label, {'color'    : '%s'%y_colour,  'fontsize'   : y_size })

    # If true, creates a legend
    if legend_on == True:
        if print_on == True:      
            print("Producing legend")

        leg = plt.legend(handles=HANDLES,labels=DATALABELS)

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
        plt.tick_params(which="both",direction="in",top=True,right=True)
    if print_on == True:
        print("Saving figure to %s"%figure_name)
    
    if save_fig == True:
        plt.savefig(figure_name)

        if print_on == True:
            print("Figure plotted! SCIENCE!")
        
    if display_figure == True:
        plt.show()
    else:
        plt.close()

# ============================================================================
#                           PARAMETERS
# ============================================================================
# FIGURE AND AXIS SIZES ======================================================
mpl.rcParams['xtick.labelsize'] = 14        # Size of x-tick labels
mpl.rcParams['ytick.labelsize'] = 14        # Size of y-tick labels
mpl.rcParams['axes.labelsize'] = 14         # Size of axes labels
mpl.rcParams['image.origin'] = 'lower'      # 

mpl.rcParams['axes.linewidth'] = 2.5        # Size of axes line width
mpl.rcParams['xtick.major.size'] = 14       # Size of major x-ticks 
mpl.rcParams['xtick.minor.size'] = 6        # Size of minor x-ticks
mpl.rcParams['xtick.major.width'] = 2.5     # Width of major x-ticks  
mpl.rcParams['xtick.minor.width'] = 1.5     # Width of minor x-ticks
mpl.rcParams['ytick.major.size'] = 14       # Size of major y-ticks
mpl.rcParams['ytick.minor.size'] = 6        # Size of minor y-ticks
mpl.rcParams['ytick.major.width'] = 2.5     # Width of major y-ticks
mpl.rcParams['ytick.minor.width'] = 1.5     # Width of minor y-ticks

# Sets font style and size
mpl.rcParams.update({'font.size': 18, 'font.weight': 'bold'})

# FIGURE PARAMETERS ==========================================================
share_plot = True               # If true, plot all variables onto same axis
                                # if false, outputs each plot in seperate
                                # figures
figure_name = 'Test.pdf'        # Defines default filename

figsize = (10,7)                # Defines figure size

display_figure = True           # Opens figure if true

save_fig = True                 # Saves the figure to file if true

print_on = True                 # Includes print statements in output if true

# POINTS & FILL PARAMETERS ===================================================
POINTSTYLES = ['*','o']         # Point style
SIZES = [1,2]                   # Point sizes
COLOURS = ['0.3','b']           # Point fill colours
EDGECOLOURS = ['face','g']      # Point edge colours

# AXIS PARAMETERS ============================================================
x_min = 0.0                 # Minimum x-axis value to plot
x_max = 10.0                # Maximum x-axis value to plot
y_min = 0.0                 # Minimum y-axis value to plot
y_max = 10.0                # Maximum y-axis value to plot

axis_range = [x_min,x_max,y_min,y_max]

x_scale = 'linear'          # Defines x-axis scaling
y_scale = 'linear'          # Defines y-axis scaling 

auto_axis = True            # Auto define dimesions of axises

both_axis = True            # Places ticks on top and right axes too

# AXIS LABELS ================================================================
x_label = "x-label"        # x-axis label
y_label = "y-label"        # y-axis label
x_size = 16                # Size of x-axis label
y_size = 16                # Size of y-axis label
x_colour = "k"             # Colour of x-axis label
y_colour = "k"             # Colour of y-axis label

minorticks = True          # Places minor ticks on axis if True
grid = False               # Places grid on plot if True

# LEGEND PARAMETERS ==========================================================
legend_on = True           # Creates legend if true

frame_border = False       # Puts border around the legend
frame = True               # Sets frame to visible
frame_colour = "k"         # Sets frame colour

# END PROGRAM ================================================================