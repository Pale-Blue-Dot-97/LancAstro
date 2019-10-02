# =====================================================================================================================
# ================================================= HISTOMAKE_TEST ====================================================
# =====================================================================================================================
#
# Example script for the use of Histomake and Plot2D modules of the LancAstro.py project.
# Adapted from code by Joshua Butterworth.
# Uses SC4K data processed by Sergio Santos and Joshua Butterworth
#
# Author: Harry Baker
# Date: 02.10.2019
# Version: 1.0
#
# =====================================================================================================================
#                                                       IMPORTS
# =====================================================================================================================

import numpy as np              # For maths functions
from astropy.io import fits     # For FITS file handling
import HistoMake as hm          # Imports HistoMake for binning data
import Plot2D as laplt          # Import Plot2D for plotting the histograms

# =====================================================================================================================
#                                               VARIABLES AND DATA
# =====================================================================================================================
# NOT IMPORTANT FOR UNDERSTANDING THE WORKINGS OF LANCASTRO
NBs = [392, 501, 711, 816]
IAs = [427, 464, 484, 505, 527, 574, 624, 679, 709, 738, 767, 827]
Zs = [2.22, 2.5, 3.1, 3.9, 4.7, 5.4]
zErrors = [0.02, 0.1, 0.4, 0.3, 0.2, 0.5]

uniID = NBs + IAs

IA_VOLUMES = [4.0e6, 4.2e6, 4.3e6, 4.3e6, 4.5e6, 4.9e6, 5.2e6, 5.5e6, 5.1e6, 5.1e6, 5.5e6, 4.9e6]
NB_VOLUMES = [0.6e6, 0.9e6, 1.2e6, 1.8e6]

zMass = []
Phis = []
Masxs = []
Comp = 1

# LOAD DATA FROM FITS FILE ============================================================================================
hdul = fits.open('DATA3.fits')
data = hdul[1].data
ID = np.array(data['ID_SC4K'])
BAND_ARRAY = [int(i[7:10]) for i in ID]  # ID=SC4K_IA427_1234 becomes BAND=IA427
uniqueID = np.unique(BAND_ARRAY)
Comments = np.array(data['Comments'])
Res = np.array(data['Response'])
AllMass = np.array(data['Mstar_1'])
Redshift = np.array(data['LAE_Redshift'])
Xray = np.array(data['XRAY_DET'])
Radio = np.array(data['RADIO_DETECTED'])
uniqueZ = np.unique(Redshift)

# =====================================================================================================================
#                                                       MAIN
# =====================================================================================================================

for j, k in zip(Zs, zErrors):
    MASK = (Comments != [' B ']) * (Comments != [' DB']) * (Redshift <= j + k) * (Redshift >= j - k) * \
           (Radio == False) * (Xray == False)
    New_mass = AllMass[MASK]
    zMass.append(New_mass)

# This will hold all the labels for each plot
labels = []

for y, x in zip(Zs, zErrors):
    z = 'z=' + str(y) + '$\pm$' + str(x)
    labels.append(z)

# Holds each set of data to be binned
DATASET = []

# Cleans data first
for M in zMass:
    CleanM = [x for x in M if str(x) != 'nan']
    DATASET.append(CleanM)

# 2D arrays that hold the x and y for each histogram plot
x_set = []
y_set = []

nbins = 10  # Number of bins to divide data into

# Uses HistoMake to bin the data and return the x and y for each histogram ready to be plotted
for i in range(len(DATASET)):
    x, y, bin_min, bin_max = hm.output_bins_to_plot(DATASET[i], n=nbins)
    x_set.append(x)
    y_set.append(y)

# Defines the colour for each plot
COLOURS = ["b", "g", "r", "c", "m", "y"]

# Uses Plot2D to plot all the x and y with the labels for the legend, the colours for the line and for the fill,
# shading the fill with alpha 0.5, defining the points to be lines, setting the x and y axis labels and the figure name
laplt.create_figure(x_set, y_set, DATALABELS=labels, COLOURS=COLOURS, FILL_COLOURS=COLOURS, POINTSTYLES=['-'],
                    x_label=r'$\log(M_\star\,/\,M_\odot)$', y_label=r'$\Phi\,/\,\mathrm{dex}^{-1}\,\mathrm{Mpc}^{-3}$',
                    SHADE=[0.5], figure_name='HistoMake_Example.pdf')
