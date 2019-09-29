import numpy as np
import ScatterPlot as sp

"""
x = np.random.rand(5,50)
y = np.random.rand(5,50)
DATALABELS = ["1","2","3","4","5"]
COLOURS = ["b","g","k","c","y"]
"""
#PATH = ["Users\hjbak\Documents\University\XGAL Internship 2019\LancAstro"]
DATANAME = ['DATA3.fits']
X_COLUMNS = [['LAE_Redshift']]
Y_COLUMNS = [['BB_mag']]

DATALABELS = ['Mass']
COLOURS = ["b"]

x = []
y = []

print("LOADING RA")
x.append(sp.data_load(DATANAME,RA))

print("\nLOADING DEC")
y.append(sp.data_load(DATANAME,DEC))

sp.create_figure(x,y,DATALABELS,COLOURS)