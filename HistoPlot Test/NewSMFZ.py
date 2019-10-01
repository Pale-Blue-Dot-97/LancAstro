import numpy as np
from astropy.io import fits
import HistoMake as hm
import Plot2D as laplt

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

for j, k in zip(Zs, zErrors):
    MASK = (Comments != [' B ']) * (Comments != [' DB']) * (Redshift <= j + k) * (Redshift >= j - k) * (
                Radio == False) * (Xray == False)
    New_mass = AllMass[MASK]
    zMass.append(New_mass)

labels = []
for y, x in zip(Zs, zErrors):
    z = 'z=' + str(y) + '$\pm$' + str(x)
    labels.append(z)

DATASET = []

for M in zMass:
    CleanM = [x for x in M if str(x) != 'nan']
    DATASET.append(CleanM)

COLOURS = ["b", "g", "r", "c", "m", "y"]

x_set = []
y_set = []

nbins = 10  # Number of bins to divide data into

for i in range(len(DATASET)):
    x, y, bin_min, bin_max = hm.output_bins_to_plot(DATASET[i], n=nbins)
    x_set.append(x)
    y_set.append(y)

laplt.create_figure(x_set, y_set, DATALABELS=labels, COLOURS=COLOURS, FILL_COLOURS=COLOURS, POINTSTYLES=['-'],
                    x_label=r'$\log(M_\star\,/\,M_\odot)$', y_label=r'$\Phi\,/\,\mathrm{dex}^{-1}\,\mathrm{Mpc}^{-3}$')
