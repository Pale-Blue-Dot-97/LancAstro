import matplotlib.pyplot as plt
import math
import pandas as pd
import numpy as np
from astropy.io import fits
#import HistoPlot as hp
import ScatterPlot as sp

NBs=[392, 501, 711, 816]
IAs=[427, 464, 484, 505, 527, 574, 624, 679, 709, 738, 767, 827]
Zs=[2.22,2.5,3.1, 3.9, 4.7, 5.4]
zErrors=[0.02,0.1,0.4,0.3,0.2,0.5]

uniID=NBs + IAs

IA_VOLUMES=[4.0e6,4.2e6,4.3e6,4.3e6,4.5e6,4.9e6,5.2e6,5.5e6,5.1e6,5.1e6,5.5e6,4.9e6]
NB_VOLUMES=[0.6e6,0.9e6,1.2e6,1.8e6]

zMass=[]
Phis=[]
Masxs=[]
Comp=1

hdul =fits.open('DATA3.fits')
data = hdul[1].data
ID = np.array(data['ID_SC4K'])
BAND_ARRAY=[int(i[7:10]) for i in ID] #ID=SC4K_IA427_1234 becomes BAND=IA427
uniqueID=np.unique(BAND_ARRAY)
Comments = np.array(data['Comments'])
Res = np.array(data['Response'])
AllMass = np.array(data['Mstar_1'])
Redshift = np.array(data['LAE_Redshift'])
Xray= np.array(data['XRAY_DET'])
Radio=np.array(data['RADIO_DETECTED'])
uniqueZ=np.unique(Redshift)
#print(uniID)
#print(uniqueZ)
for j,k in zip(Zs,zErrors):
 #print(j)
 MASK =(Comments!=[' B '])*(Comments!=[' DB'])*(Redshift<=j+k)*(Redshift>=j-k)*(Radio==False)*(Xray==False)
 New_mass=AllMass[MASK]
 zMass.append(New_mass)
 #print(New_mass)

#labels=[]
#for y,x in zip(Zs,zErrors):
# z='z=' + str(y) + '$\pm$' + str(x)
# labels.append(z)

DATASET = []
#COLOURS = ["b","g","r","c","m","y"]
COLOURS = ['tab20b',"green"]

for M in zMass:
 #M     = np.array(data['Mstar_1'])       #Read stellar masses in logMsun
 CleanM = [x for x in M if str(x) != 'nan']
 DATASET.append(CleanM)

 #print("Printing dataset")
 #print (CleanM)
 #hp.create_figure(DATASET,labels,COLOURS)

 nbins = 10                            #Number of bins to divide data into
 V     = 1e5                             #Survey volume in Mpc3
 Phi,edg = np.histogram(CleanM,bins=nbins) #Unnormalized histogram and bin edges
 #print(edg)
 dM    = edg[1] - edg[0]                 #Bin size
 Masx   = edg[0:-1] + dM/2.               #Mass axis
 Phi   = Comp*Phi / V / dM                    #Normalize to volume and bin size
 #print(Phi)
 #MASK2=(Phi !=0.)
 #Masx=Masx[MASK2]
 #Phi= Phi[MASK2]
 #CleanPhi = [x for x in M if str(x) != 'nan']
 Phis.append(Phi)
 Masxs.append(Masx)
 #Create the Respective labels
labels=[]
for y,x in zip(Zs,zErrors):
 z='z=' + str(y) + '$\pm$' + str(x)
 labels.append(z)

print(len(DATASET))
print(len(Phis))

sp.create_figure(Masxs,Phis,labels)
""" 
for i in np.arange(0, len(labels),1):
 Ms=Masxs[i]
 Ps=Phis[i]
 Ls=labels[i]
 plt.plot(Ms,Ps, marker='o',label=Ls)
plt.yscale('log')
plt.xlabel(r'$\log(M_\star\,/\,M_\odot)$')
plt.ylabel(r'$\Phi\,/\,\mathrm{dex}^{-1}\,\mathrm{Mpc}^{-3}$')
plt.tight_layout()
plt.legend()
plt.show()
"""