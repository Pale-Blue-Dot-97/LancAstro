############################## HISTOPLOT ##############################
#
# Script to plot histograms of data stored in FITS files, ascii tables
#
# Author: Harry Baker
# Date: 10.07.2019
# Version: 1.3
# ==VERSION 1==
#
# Part of the LancAstro.py project for XGAL
#
# Edited from HISTO_Environment.py by Dr David Sobral 2008
#
# TASKS:
# -Add ASCII read capability
# -Improve default settings
# -Improve random colour assignment
# -Improve methods to take more optional arguements
# -Fully comment
#
version = "1.3"
############################### IMPORTS ################################

import os, sys
import matplotlib.pyplot as plt
from astropy.io import fits as pyfits
import numpy as np
from mpltools import color  
import random
import Dependencies.ascii_read as ar
import Dependencies.progress_bar as pb

################################# DEFINE DATA ##########################

DATANAME = ['DATA3.fits']            # Data filenames
COLUMNNAME = [['Mstar_1']]           # Name of FITS Column name    
DATALABEL = ["Mass"]                 # Label of each variable for legend
PATH = ["",""]                       # Path to data files

################################# FILENAME ##############################
figure_name = "Josh_Test.png"    # Name of plot to write to

########################### PARAMETERS ###################################

# LINE & FILL PARAMETERS =============================================
LINES = ['-','-']           # Line style
LINW = [1,2]                # Line width
COLOURS = ['0.3','b']       # Line colours
SHADE = [0.6,0.5]

# BIN PARAMETERS ===================================================
WID = [1.0,1.0]             # Width of bins
MIN = [7.0,7.0]             # Minimum value
MAX = [15.0,10.0]           # Maximum value
n = [150,10]                # Number of bins

use_width = False          # Defines if WID is used (True) to construct bins or n is used (False)
auto_bins = True           # Auto create MIN, MAX, WID 

# AXIS DIMENSIONS ==================================================
x_min = 0.0                 # Minimum x-axis value to plot
x_max = 15.0                # Maximum x-axis value to plot
y_min = 0.0                 # Minimum y-axis value to plot
y_max = 20.0                # Maximum y-axis value to plot

auto_axis = True            # Auto define dimesions of axises

# AXIS LABELS ======================================================
x_label = r"Log(M/M$\odot$)"    # x-axis label
y_label = "Counts"              # y-axis label
x_size = 16                     # Size of x-axis label
y_size = 16                     # Size of y-axis label
x_colour = "k"                  # Colour of x-axis label
y_colour = "k"                  # Colour of y-axis label

minorticks = True          # Places minor ticks on axis if True
grid = False               # Places grid on plot if True

# LEGEND PARAMETERS ================================================
frame_border = False       # Puts border around the legend
frame = True               # Sets frame to visible
frame_colour = "k"         # Sets frame colour

############################ DATA LOADING #################################
DATASET = []

def data_load(PATH, DATANAME, COLUMNNAME):

    pb.printProgressBar(0,len(DATANAME), prefix = 'Loading Data:', suffix = 'Complete')
    for i in range(len(DATANAME)):
        print("Opening Data file: %s%s"%(PATH[i],DATANAME[i]))
        try:
            FILE = pyfits.open(PATH[i] + DATANAME[i])
            pb.printProgressBar(i,len(DATANAME), prefix = 'Loading Data:', suffix = 'Complete')
        except:
            print("FILE COULD NOT BE OPENED")

        TABLE = FILE[1].data
        for j in COLUMNNAME[i]:
            print("Loading column %s"%j)
            VAR = TABLE.field("%s"%j)
            CleanVAR = [x for x in VAR if str(x) != 'nan']
            DATASET.append(CleanVAR)

        print("Data set loaded ready for plotting")
    return DATASET

def convert_to_array(var):
    a = []
    if type(var) == 'float' or 'double':
        a.append(var)
        return a

    else:
        return var

###################### FIND MIN/MAX OF DATA ################################
def _get_outer_edges(a):
    """
    Determine the outer bin edges to use, from either the data or the range
    argument
    ########
    Edited from Numpy source code:
    https://github.com/numpy/numpy/blob/master/numpy/lib/histograms.py
    """
    #print("Determining edges of data")
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

################################# CONSTRUCT BINS ##########################################

def make_bins(a,MIN,MAX,WID,n):
	# Calculates the outer edges of the bins
    if auto_bins == True:
       MIN, MAX = _get_outer_edges(a)

    #print("Determining width of bins from requested number of bins")
    if use_width == False:
        WID = (MAX - MIN)/float(n)

    #print("Constructing bins")
	# Constructs the centre of the bins from the defined mid, max and width
    bins = np.arange(MIN+WID/2.0,MAX-WID/2.0,WID)
 
	# Creates an array of the same length as bins to hold the counts per bin
    counts = []
    counts = bins - bins

	# Loops through all data points, for each one cycling through the bins to
	# see if it resides within it. If so, it is added to counts for that bin
    for j in range(len(a)):
        for k in range(len(bins)):
                if a[j] < bins[k] + WID/2.0:
                    if a[j] > bins[k] - WID/2.0:
                        counts[k] = counts[k] + 1.0

	# Normalises the data
    NORM = sum(counts)

	# Constructs the extra points either side of the mid-point of the bin to
	# plot
    x = [MIN] # Array of x-axis plot points (bin edges)
    y = [y_min] # Array of y-axis plot points (count values for each bin)
    for j in range(len(bins)):

        # Left edge
        x.append(bins[j] - WID/2.0)
        y.append(counts[j])

		# Central value
        x.append(bins[j])
        y.append(counts[j])

		# Right edge
        x.append(bins[j] + WID/2.0)
        y.append(counts[j])

    x.append(MAX)
    y.append(y_min)

    return x, y, NORM, MIN, MAX

#########################################################################

"""
Checks every parameter list is the same length
i.e if user only defines min and max for 1 variable
Corrects these to all have same length  
"""
def default_lengths(param, DATALABEL):
    if len(param) != len(DATALABEL):
        print("WARNING: Not every variable has defined range limits!")
        print("Setting all variables to default of 1st entry of MIN, MAX, WID, n \n")
        for i in range(len(param)):
            param[i] = param[0]

        for i in range(len(DATALABEL) - len(param)):
            param.append(param[0])

    return param

def colour_cycle(n_var,ax,COLOURS):
    if len(COLOURS) != n_var:
        for i in range(n_var):
            r = lambda: random.randint(0,255)
            random_colour = '#%02X%02X%02X' % (r(),r(),r())
            if i < len(COLOURS):
                COLOURS[i] = random_colour
            else:
                COLOURS.append(random_colour)
        #COLOURS = color.cycle_cmap(n_var, cmap='%s'%COLOURS[0])
        #COLOURS = plt.get_cmap("%s"%random_colour)
        #cmap = plt.get_cmap("%s"%COLOURS[0])
        #print(cmap)

#################################################################################################

def plot(DATASET,DATALABEL,COLOURS):
    PLOTS = []
    HANDLES = []
    LEGEND = []
    all_counts = []

    global MIN, MAX, WID, n, LINES, LINW, SHADE
    MIN = default_lengths(MIN,DATALABEL)
    MAX = default_lengths(MAX,DATALABEL)
    WID = default_lengths(WID,DATALABEL)
    n = default_lengths(n,DATALABEL)
    #COLOURS = default_lengths(COLOURS,DATALABEL)
    LINES = default_lengths(LINES,DATALABEL)
    LINW = default_lengths(LINW,DATALABEL)
    SHADE = default_lengths(SHADE,DATALABEL)

    colour_cycle(len(DATALABEL),PLOTS,COLOURS)

    pb.printProgressBar(0,len(DATASET),prefix='BEGINNING PLOTTING',suffix='COMPLETE',length=40)

    for i in range(len(DATASET)):
        #print("Plotting Variable %d"%(i+1))

        pb.printProgressBar(i,len(DATASET),prefix='PLOTTING VARIABLE %d'%(i+1),suffix='COMPLETE',length=40)

        x, y, NORM, MIN[i], MAX[i] = make_bins(DATASET[i],MIN[i],MAX[i],WID[i],n[i])
        
        all_counts.append(y)
        
        LINE = plt.plot(np.array(x),np.array(y), linestyle=LINES[i],color=COLOURS[i],linewidth=int(LINW[i]), label = None) 
        FILL = plt.fill(np.array(x),np.array(y), alpha=SHADE[i],color=COLOURS[i],label=DATALABEL[i])
        HANDLES.append(FILL[0])
    
    pb.printProgressBar(len(DATASET),len(DATASET),prefix='FINISHED PLOTTING!',length=40)    
    
    return HANDLES, all_counts, MIN, MAX

#####################################################################

def determine_axis(MIN,MAX,all_counts):
    #print("Determing dimensions of axis")
    x_min = np.min(MIN)
    x_max = np.max(MAX)

    max_counts = []
    for i in all_counts:
        max_counts.append(np.max(i))

    y_max = 1.1 * np.max(max_counts)

    return x_min, x_max, y_max

#############################################################################

def create_figure(DATASET,DATALABEL,COLOURS=["b"],figure_name="Fig1.png"):
    print("WELCOME TO HistoPlot %s"%version)
    print("PART OF THE LancAstro PACKAGE \n")

    # Calls plot method to plot each variable from the dataset
    # plot() calls methods to construct bins    
    HANDLES, all_counts, MIN, MAX = plot(DATASET,DATALABEL,COLOURS)

    # If set to, automatically sets the min and max of the axis
    # from the range of all the data
    if auto_axis == True:
        x_min, x_max, y_max = determine_axis(MIN,MAX,all_counts)

    # Places grid if set to do so
    plt.grid(grid)

    # Adds minor ticks to axis if set to do so
    if minorticks == True:
        plt.minorticks_on()

    # Creates figure axis from defined parameters
    plt.axis([x_min,x_max,y_min,y_max])

    # Creates x and y axis labels from settings 
    plt.xlabel(r'%s'%x_label, {'color'    : "%s"%x_colour,  'fontsize'   : x_size })
    plt.ylabel(r'%s'%y_label, {'color'    : '%s'%y_colour,  'fontsize'   : y_size })

    print("Producing legend")
    legend = plt.legend(HANDLES,DATALABEL)

    # Sets legend frame to be transparent if set to do so
    if frame == False:
        legend.get_frame().set_facecolor('none')

    # Removes legend border if set to False
    if frame_border == False:
        legend.get_frame().set_linewidth(0.0)

    # Sets frame border to defiend
    if frame_border == True:
        legend.get_frame().set_edgecolor(frame_colour)
    
    print("Saving figure to %s"%figure_name)
    plt.savefig(figure_name)

    print("Figure plotted! SCIENCE!")
    plt.show()

#================================ TESTING ====================================
# Loads the data into Python
#DATASET = data_load(PATH, DATANAME, COLUMNNAME)

# Calls in methods to create bins, plot. display and save figure
#create_figure(DATASET, DATALABEL, COLOURS)

# END PROGRAM