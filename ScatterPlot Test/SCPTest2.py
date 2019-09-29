import os
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits as pyfits
import Dependencies.ascii_read as ar
import Dependencies.progress_bar as pb
import ScatterPlot as sp 
import Plot2D as laplt
band_nums = ["427","464","484","505","527","574","624","679","709","711","738","767","827"]
path = ["Plot_data\\"]

filenames = []

for i in band_nums:
	filenames.append("Plot_data_Band%s.fits"%i)

columnnames = ["Log_Mass_Center","Log_Phi_no_completeness","Log_Phi_Error_no_completeness_minus","Log_Phi_Error_no_completeness_plus"]

columns = []

for i in range(len(filenames)):
	columns.append(columnnames)

data = laplt.data_load(filenames,columns,path)

x = []
y = []
y_err_lower = []
y_err_upper = []

for i in range(len(filenames)):
	x.append(data[4*i])
	y.append(data[4*i + 1])
	y_err_lower.append(data[4*i + 2])
	y_err_upper.append(data[4*i + 3])

y_err = []

for i in range(len(y_err_lower)):
	y_err.append([y_err_lower[i],y_err_upper[i]])

pointstyle = ['o']

laplt.create_figure(x,y,y_error=y_err,DATALABELS=band_nums,POINTSTYLES=pointstyle)
