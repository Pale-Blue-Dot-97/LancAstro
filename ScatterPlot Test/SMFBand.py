import matplotlib.pyplot as plt
plt.style.use('seaborn-whitegrid')
import matplotlib.cm as cm
import math
import pandas as pd
import numpy as np
from astropy.io import fits
import ScatterPlot as sp 

NBs=np.array([392, 501, 711, 816]) # NB filter band numbers
IAs=np.array([427, 464, 484, 505, 527, 574, 624, 679, 709, 711, 738, 767, 827]) # IA filter band numbers
Zs=[2.22,2.5,3.1, 3.9, 4.7, 5.4] # All Redshift values
zErrors=[0.02,0.1,0.4,0.3,0.2,0.5] # Redshift error
newZs=[2.5,3.1, 3.9, 4.7, 5.4] # reduced Redshift values, neglecting 2.22 as it causes an error when the NB filters are not used
newzErrors=[0.1,0.4,0.3,0.2,0.5] # Redshift range
uniID=np.concatenate((NBs,IAs),axis=None) #Create an array of both NB and IA band values
IA_VOLUMES=np.array([4.0e6,4.2e6,4.3e6,4.3e6,4.5e6,4.9e6,5.2e6,5.5e6,5.1e6,1.2e6,5.1e6,5.5e6,4.9e6]) # Surveyed volume in IA filters
NB_VOLUMES=np.array([0.6e6,0.9e6,1.2e6,1.8e6]) # Surveyed volume in NB volume
uniVol=np.concatenate((NB_VOLUMES,IA_VOLUMES),axis=None) # #Create an array of both NB and IA surveyed volume values




Comp=1
hdul =fits.open('DATA3.fits')
data = hdul[1].data
ID = np.array(data['ID_SC4K'])
BAND_ARRAY=[int(i[7:10]) for i in ID] #ID=SC4K_IA427_1234 becomes BAND=427
BAND_ARRAY=np.array(BAND_ARRAY)
uniqueID=np.unique(BAND_ARRAY)
Comments = np.array(data['Comments'])
Res = np.array(data['Response'])
AllMass = np.array(data['Mstar_1'])
Redshift = np.array(data['LAE_Redshift'])
Xray= np.array(data['XRAY_DET'])
Radio=np.array(data['RADIO_DETECTED'])
Flux=np.array(data['Flux'])

Flux=[math.log(x,10) for x in Flux]
Flux=np.array(Flux)
Phinocos=[]
uniqueZ=np.unique(Redshift)
zMass=[] # Initialise the Array which will contain filtered mass values. 
Phis=[]
PhiErrors=[]
CoPhiErrors=[]
Masxs=[]
Masx=[]
zID=[]
zF=[]
edg=[]
#print(uniID)
#print(uniqueZ)
######################################################################################################
def make_bins(a,N,ID,F):
    
    #Calculate max and in values of the array
    MAX=np.amax(a)
    MIN=np.amin(a)
    # Calculates the width of the bins
    WID=(MAX-MIN)/(N)
    
    
    CENT_MIN=MIN+WID/2.
    CENT_MAX=MAX-WID/2.
    #print("Constructing bins")
    # Constructs the centre of the bins from the defined mid, max and width
    bins = np.arange(CENT_MIN,CENT_MAX+WID,WID)
    MASK=(bins<=MAX)*(bins>=MIN)
    
    bins=bins[MASK]
    
    edg = np.arange(MIN,MAX+WID,WID)
    #MASK2=(len(edg)<=len(bins)+1)
    #edg=edg[MASK2]
    if len(edg)>N+1:
       edg=edg[:-1]
    uID=np.unique(ID)
    countsFilt=[]
    FiltsVol=[]
    
    for u in uID:
        #print(u)
        MASK3=(ID==u)
        MASK4= (IAs==u)
        FiltID=ID[MASK3]
        FiltVol=IA_VOLUMES[MASK4]
        countFilt=len(FiltID)
        countsFilt.append(countFilt)
        FiltsVol.append(FiltVol[0])      
    
    Vol=np.sum([x*y for x,y in zip(countsFilt,FiltsVol)])/sum(countsFilt)
    #print(Vol)
    # Creates an array of the same length as bins to hold the counts per bin
    counts = []
    counts = bins - bins
    countsComp = []
    countsComp = bins - bins
    # Loops through all data points, for each one cycling through the bins to
    # see if it resides within it. If so, it is added to counts for that bin
    for j in range(len(a)):
        for k in range(len(bins)):
                if a[j] <= bins[k] + WID/2.0:
                    if a[j] >= bins[k] - WID/2.0:
                        iD=ID[j]
                        flu=F[j]
                        hdul2 =fits.open('COMPLETENESS_functions/COMPLETENESS_IA%d.fits'%iD)
                        data2 = hdul2[1].data
                        flux = np.array(data2['lineflux'])
                        comp = np.array(data2['completeness'])
                        comp= np.interp(flu, flux, comp)
                        comp=comp/100
                        countsComp[k] = countsComp[k] + 1/comp
                        counts[k] = counts[k] + 1
    # Normalises the data
    NORM = sum(counts)

    # Constructs the extra points either side of the mid-point of the bin to
    # plot
    x = [] # Array of x-axis plot points (bin edges)
    co = [] # Array of y-axis plot points (count values for each bin)
    z = [] # Array of central x-axis points
    cocomp = []
    for j in range(len(bins)):
        # Central value
        co.append(counts[j])
        cocomp.append(countsComp[j])
    return bins, co, cocomp, edg, WID, NORM, MIN, MAX, Vol
#########################################################################################



for j in uniID:
    #print(j)
    MASK =(Comments!=[' B '])*(Comments!=[' DB'])*(Radio==False)*(Xray==False)*(BAND_ARRAY==j)*(BAND_ARRAY!=392)*(BAND_ARRAY!=501)*(BAND_ARRAY!=816)
    New_mass=AllMass[MASK]
    New_ID=BAND_ARRAY[MASK]
    #print(New_ID)
    New_flux=Flux[MASK]
    if len(New_ID)>0:
		
        zMass.append(New_mass)
        zID.append(New_ID)
        zF.append(New_flux)
#print(zID)


for M,I,F in zip(zMass,zID,zF):
    Phiboth=[]
    Mask=(str(M)!='nan')   
    #M     = np.array(data['Mstar_1'])       #Read stellar masses in logMsun
    #CleanM = [x for x in M if str(x) != 'nan']
    CleanM=M[Mask]
    CleanI=I[Mask]
    CleanF=F[Mask]
    useI=CleanI[0]
    useF=CleanF[0]
    useM=CleanM[0]
    #print(useM)
    nbins = 10                         #Number of bins to divide data into
    Masx, Phinoco,Phi,edg, Wid, Norm, Min, Max, Vol=make_bins(useM, nbins,useI,useF) 
    Phi=np.array(Phi)
    #print(Vol)
    Phinoco=np.array(Phinoco)
    #print(Phinoco)
    PhiError=np.array([math.sqrt(x) for x in Phinoco])
    PhiError= PhiError / Vol / Wid  
    Phi   = Phi / Vol / Wid   #Normalize to volume and bin size
    Phinoco   = Phinoco / Vol / Wid
    

    #print(Phi)
    #Cofa= [x/y for x,y in zip(Phi,Phinoco)]
    #PhiError= np.array([x*y for x,y in zip(PhiError,Cofa)])

 
    
    

    ##########################################
    #Apply MASK2 to output arrays to eliminate 0 values.
    MASK2=(Phi > 0.)
    PhiError=PhiError[MASK2]
    Masx=Masx[MASK2]
    Phi= Phi[MASK2]
    Phinoco = Phinoco[MASK2]
    Cofa= [x/y for x,y in zip(Phi,Phinoco)]
    
    CoPhiError= np.array([x*y for x,y in zip(PhiError,Cofa)])
    #print([x/y for x,y in zip(Phi,Phinoco)])
    Phis.append(Phi)
    PhiErrors.append(PhiError)
    CoPhiErrors.append(CoPhiError)
    Phinocos.append(Phinoco)
    Masxs.append(Masx)
 #Create the Respective labels
labels=[]
labels2=[]
#print(PhiErrors)
for x in uniID:
    z=str(x)
    a=str(x) + ' (Without Completeness)'
    labels.append(z)
    labels2.append(a)


colors= cm.rainbow(np.linspace(0,1, len(Phis)))
colors2= cm.rainbow(np.linspace(1,1, len(Phis)))

size = [5.0]

#for i in range(len(Masxs)):
#    size.append(size[0])

sp.create_figure(Masxs,Phinocos,y_error=PhiErrors,DATALABELS=labels,COLOURS=colors,SIZES=size)

"""
for i,c,c2,B in zip(np.arange(0, len(labels),1),colors,colors2,uniID):
    Ms=Masxs[i]
    Ps=Phis[i]
    PErrs=PhiErrors[i] 
    coPErrs=CoPhiErrors[i]            
    Ls=labels[i]
    Ls2=labels2[i]
    Phinc=Phinocos[i]
    check=np.array([x-y for x,y in zip(Ps,coPErrs)])
    #MASK4=(check==0.)
    
    #print(check)
    #print(Ps)
    #print(coPErrs)
    #print(Ps-coPErrs)
    logP=np.array([math.log(x,10) for x in Ps])
    for i,j in zip(check,np.arange(0,len(check),1)):
        if i <= 0.:
            check[j]=10**(logP[j]-1)
            
    #print(check)
    logPhinc=[math.log(x,10) for x in Phinc]
    cologPp=[math.log(x+y,10) for x,y in zip(Ps,coPErrs)]
    cologPm=[math.log(x,10) for x in check]
    minco=[x-y for x,y in zip(logP,cologPm)]
    Plusco=[x-y for x,y in zip(cologPp,logP)]
    coasy_err=[minco,Plusco]
    #print(coasy_err)
    check=np.array([x-y for x,y in zip(Phinc,PErrs)])
    MASK4=(check==0.)
    
    logPp=[math.log(x+y,10) for x,y in zip(Phinc,PErrs)]
    for i,j in zip(check,np.arange(0,len(check),1)):
        if i == 0.:
            check[j]=10**(logPhinc[j]-1)
    logPm=[math.log(x,10) for x in check]
    minnoco=[x-y for x,y in zip(logPhinc,logPm)]
    Plusnoco=[x-y for x,y in zip(logPp,logPhinc)]
    asy_err=[minnoco,Plusnoco]
    #print(asy_err)
    #print(coasy_err)
    #plt.plot(Ms,logP,'o',label=Ls,color=c)
    #plt.plot(Ms,logP,'o',label=Ls,color=c)
    plt.errorbar(Ms,logP,yerr=coasy_err,fmt='o',color=c,label=Ls,alpha=0.5)
    plt.errorbar(Ms,logPhinc,yerr=asy_err,fmt='o',color=c2,label=Ls2,alpha=0.5)
    #plt.errorbar(Ms,Phinc,fmt='o',color=c2,label=Ls)
    #plt.errorbar(Ms,Ps,yerr=PErrs,fmt='o',color=c,label=Ls)
    #plt.fill_between(Ms,Ps-PErrs,Ps+PErrs)
    #plt.yscale('log', nonposy='clip')
    plt.xlabel(r'$\log(M_\star\,/\,M_\odot)$')
    plt.ylabel(r'$\log(\Phi\,/\,\mathrm{dex}^{-1}\,\mathrm{Mpc}^{-3})$')
    #plt.ylabel(r'$\Phi\,/\,\mathrm{dex}^{-1}\,\mathrm{Mpc}^{-3}$')
    #plt.ylim(1.e-10,1.e-5)
    #plt.autoscale(enable=True,axis='y',tight=True)
    plt.legend()
    plt.tight_layout()
    plt.savefig('FiltSMFs/SMF_Redshift%s.png'%B)
    plt.clf()
"""