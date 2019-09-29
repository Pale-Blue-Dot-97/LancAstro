
# ============================================================================
#                             IMPORTS   
# ============================================================================

import os, sys
import matplotlib as mpl
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import numpy as np 
import Dependencies.progress_bar as pb
import Dependencies.DataLoad as dl
import Plot3D as laplt

#=============================================================================
filename = ['JOSH4_SC4K_LAEs_catalogue_Reduced_sorted.fits']
columnnames = [['RA','Dec','LAE_Redshift']]

data = dl.data_load(filename,columnnames)

RA = data[0]
Dec = data[1]
z = data[2]

cmap = cm.get_cmap('jet')
colour = 'colourbar'

laplt.create_visualisation(RA,Dec,z=z,DATALABELS=columnnames,POINTSTYLES='o',
					figsize=(12,8),figure_name='SC4K_Data_Cube',x_label='RA',
					y_label='Dec',z_label='Redshift',CMAP=[cmap],COLOURS=[colour],
					no_frames=50,elev=90)